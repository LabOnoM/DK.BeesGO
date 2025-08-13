---
title: "20250729:fenA processx_convert & theme"
draft: false
tags:
  - Dry
  - Linux
  - ShinyApp
  - Functional-Enrichment-Analysis
---
# Today's task

1. Make fenA fit with our [loc](https://github.com/LabOnoM/ShinyAppsManager) platform
2. Convert the main loop of fenA ShinyApp into Subprocesses of R session by using [processx](https://github.com/r-lib/processx)
3. Change the theme of fenA ShinyApp into the uniformed format
4. If possible, draw software scheme for the main logic in fenA ShinyApp

# 1. Add code blocks for [loc](https://github.com/LabOnoM/ShinyAppsManager)

## Basic functions
1. Handle crash (inside **header** of navbarPage)
2. Handle shiny theme selection
3. Handle live console
4. Email notification (inside **header** of navbarPage)
5. Live Console (inside **footer** of navbarPage)
6. Communications with [loc](https://github.com/LabOnoM/ShinyAppsManager) based on [jwt](https://en.wikipedia.org/wiki/JSON_Web_Token)(inside **header** of navbarPage)

Helper functions (global.R)
```R
# Get Info ----
print_to_console <- function(msg) {
  js_line <- gsub("\\\\", "\\\\\\\\", msg)
  js_line <- gsub("\"", "\\\\\"", js_line)
  shinyjs::runjs(sprintf(
    "if (document.getElementById('console')) {
       let el = document.getElementById('console');
       let lines = el.innerText.split('\\n');
       if (lines.length > 2000) {
         lines = lines.slice(lines.length - 1999); // keep last 1999 lines
       }
       lines.push(\"%s\");
       el.innerText = lines.join('\\n');
       el.scrollTop = el.scrollHeight;
     } else {
       console.warn('#console element not found!');
     }",
    js_line
  ))
  return(js_line)
}
# 获取当前登录用户邮箱
get_user_email <- function(token){
  headers <- add_headers(`Token` = token, `Content-Type` = "application/json")
  res <- try(GET('http://10.2.26.152:7156/sqx_fast/app/user/info', headers), silent = TRUE)
  if(!inherits(res, "try-error") && status_code(res) == 200){
    res_data <- fromJSON(content(res, "text", encoding = "UTF-8"))
    return(res_data$data$email)
  }
  return(NA)
}
if(!is.na(Sys.getenv("SHINY_SERVER_VERSION", unset = NA))){
  smtp_creds <- creds_envvar(
    user = Sys.getenv("EMAIL_USER"),
    pass_envvar = "EMAIL_PASS",
    host = "smtp.gmail.com",
    port = 587,
    use_ssl = TRUE
  )
  email <- "Email_Template.html" %>% file.path("..", "GlobalEnvironment", .) %>% readLines(., , warn = FALSE) %>% 
    paste(., collapse = "\n") %>% HTML() %>% compose_email(body = .)
  email_fail <- "Email_Template_fail.html" %>% file.path("..", "GlobalEnvironment", .) %>% readLines(., , warn = FALSE) %>% 
    paste(., collapse = "\n") %>% HTML() %>% compose_email(body = .)
  # get user home
  UsrHome <- file.path("/mnt","Public")
}else{
  UsrHome <- path_home() %>% as.character()
}
options(shiny.fullstacktrace = TRUE)
```

Example codes: Error Trigger (server.R):
```R
tryCatch(
    ., error = function(e) {
      isolate(values$error_state <- paste0("Error:\n", e$message))
      showNotification(e$message, type = "error")
      return(NULL)
    })
```

Example codes (server.R):
```R
  # Initiation ----
  px <- reactiveVal(NULL)
  # Make reactive values ----
  values<-reactiveValues(
    is_on_shiny_server = !is.na(Sys.getenv("SHINY_SERVER_VERSION", unset = NA))
  )
  para<-reactiveValues()
  
  # Fetch information if running on BSGOU loc (local allocation controller) system.
  observe({
    req(values$is_on_shiny_server)
    req(values$server_connected)
    req(is.null(session$userData$email))
    invalidateLater(10000, session)
    msg <- "Checking email"
    isolate(values$console_log <-c(values$console_log, print_to_console(msg)))
    isolate({values$console_log <- c(values$console_log, msg)})
    session$userData$email <- get_user_email(session$userData$token)
    isolate(UserEmail <- session$userData$email)
    msg <- paste0("Your email is: ", UserEmail)
    isolate(values$console_log <-c(values$console_log, print_to_console(msg)))
    isolate({values$console_log <- c(values$console_log, msg)})
    msg <- "Welcome to BSGOU!"
    isolate(values$console_log <-c(values$console_log, print_to_console(msg)))
    isolate({values$console_log <- c(values$console_log, msg)})
    isolate(values$user_email_checked <- TRUE)
  })
  output$EmailNotif_ui <- renderUI({
    req(values$user_email_checked )
    req(!is.null(session$userData$email))
    req(values$is_on_shiny_server)
    label <- paste0("Enable Email Notification to: ", session$userData$email)
    checkboxInput("EmailNotif", label = label, value = TRUE)
  })
  observe({
    msg <- "The console is working"
    isolate(values$console_log <-c(values$console_log, print_to_console(msg)))
    req(isolate(values$is_on_shiny_server))
    ## start-3嵌入代码开始，作用：进入/退出HTTP请求记录
    query <- parseQueryString(session$clientData$url_search)
    session$userData$id <- query$id
    session$userData$appName <- query$appName
    session$userData$token <- query$token
    headers <- httr::add_headers(`Token`=session$userData$token, `Content-Type`="application/json")
    connect_req = list(`appName`=session$userData$appName, `action`="connect", `id`=session$userData$id)
    connect_data <- try(
      url_execute(
        2, 'http://10.2.26.152/sqx_fast/app/workstation/shiny-connect-action',
        toJSON(connect_req, pretty = FALSE,auto_unbox = TRUE), headers), silent = TRUE
    )
    values$server_connected<-TRUE
    if (connect_data$code!=0) {session$close()}
    ## end-3嵌入代码开始，作用：进入/退出HTTP请求记录
  })
  output$conditional_head <- renderUI({
    if(isolate(values$is_on_shiny_server)){
      ## start-1嵌入代码开始，作用：异常跳转到预约系统首页
      tags$head(
        tags$script(HTML("
        $(document).on('shiny:disconnected', function(event) {
          window.location.href = 'http://10.2.26.152/';
        });
      "))
      )
      ## end-1嵌入代码结束，作用：异常跳转到预约系统首页
    }else{
      NULL
    }
  })
  ## Handle crash info ----
  # Pop up dialog to show error message
  observeEvent(values$error_state, {
    if(!is.null(values$error_state)){
      SoftwareCrashInfo(error_state = values$error_state, input = input, output = output, main_session=session)
      values$error_state<-NULL
    }
  })
  # Observe button click to close modal
  observeEvent(input$Crash_modal_ok, {
    shiny::removeModal()
  })
  # Handle shiny theme selection ----
  observe({
    req(is.null(para$UI_Init))
    updateSelectInput(
      session = session, inputId = "shinytheme-selector",
      selected = "yeti"
    )
    isolate(para$UI_Init <- "Done")
  })
  # Handle live console----
  # Clean console history
  observeEvent(input$clear_console, {
    runjs("
    var consoleEl = document.getElementById('console');
    if (consoleEl) {
      consoleEl.innerText = '';
    } else {
      console.warn('#console element not found!');
    }
  ")
    isolate({values$console_log <- NULL})
  })
  # Download console monitor log
  output$download_console <- downloadHandler(
    filename = function() {
      paste0("console_log_", Sys.Date(), ".log")
    },
    content = function(file) {
      writeLines(values$console_log, file)
    }
  )


  # Session End Function ----
  # Shutdown everything if connection closed
  session$onSessionEnded(function() {
    if(!is.na(Sys.getenv("SHINY_SERVER_VERSION"))){
      headers <- add_headers(`Token`=session$userData$token, `Content-Type`="application/json")
      connect_req_end = list(`appName`=session$userData$appName, `action`="disconnect", `id`=session$userData$id)
      connect_end_data <- try(
        url_execute(
          2, 'http://10.2.26.152/sqx_fast/app/workstation/shiny-connect-action',
          toJSON(connect_req_end, pretty = FALSE, auto_unbox = TRUE), headers
        ), silent = TRUE
      )
    }
    stopApp()
  })
```

Example codes (ui.R):
```R
        ## For BSGOU ----
        header = tagList(
          ### For lac system----
          uiOutput("conditional_head"), ## start-1嵌入代码开始，作用：异常跳转到预约系统首页
          ## Costume CSS----
          tags$style(HTML(
            '.navbar-nav > li > a, .navbar-brand {
                   padding-top: 3px !important; 
                   padding-bottom: 0px !important;
                   margin-top: 2px;
                   margin-bottom: 0px;
                   height: 25px;
                 }
             .navbar {
                min-height:26px !important;
                z-index:99999 !important;
                margin-top: 1px;
                margin-bottom: 1px;
             }'
          )),
          tags$head(
            ## Google Analytics----
            # Embed Google Analytics
            tags$script(async = NA, src = "https://www.googletagmanager.com/gtag/js?id=G-RWF876HNKW"),
            tags$script(HTML("
              window.dataLayer = window.dataLayer || [];
              function gtag(){dataLayer.push(arguments);}
              gtag('js', new Date());
              gtag('config', 'G-RWF876HNKW');
            ")),
            ## Favicons----
            tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "GREP1_Logo_32x32.png"),
            tags$link(rel = "icon", type = "image/png", sizes = "16x16", href = "GREP1_Logo_32x32.png")
          ),
          h5(HTML(
            "<strong>B</strong>ulk <strong>R</strong>NA-seq <strong>A</strong>uto-<strong>P</strong>ipeline"
          ))
        ),
        ## Live Console----
        footer = absolutePanel(
          id = "ConsoleRealTime", class = "panel panel-default",
          bottom = 50, left = 25, width = "75%", fixed = FALSE,
          draggable = TRUE, height = "100px",
          style = "opacity: 0.95; z-index: 20;", ## z-index modification
          uiOutput("EmailNotif_ui"),
          # Horizontal bar
          tags$div(
            style = "display: flex; align-items: center; gap: 10px; padding-left: 5px;",
            HTML('<button data-toggle="collapse" data-target="#GlobalSettingMainPanel">Collapse</button>'),
            tags$a("Console Monitor"),
            downloadButton("download_console", "💾 Download Console Log"),
            actionButton("clear_console", "🧹 Clear Console", icon = icon("eraser"), class = "btn-danger")
          ),
          tags$div(
            id = "GlobalSettingMainPanel",
            style = "max-height: 300px; overflow-y: auto; padding: 5px;",
            tags$div(
              id = "console",
              style = "white-space: pre-wrap; font-family: monospace; background-color: #111; color: #0f0; height: 280px; overflow-y: scroll; padding: 10px; border: 1px solid #333;"
            )
          )
        ),
```
# 2. Convert Main loop with [processx](https://github.com/r-lib/processx)

