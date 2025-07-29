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
1. Handle crash
2. Email notification
3. Communications with [loc](https://github.com/LabOnoM/ShinyAppsManager) based on [jwt](https://en.wikipedia.org/wiki/JSON_Web_Token)

Example codes:
```R
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
```

# 2. Convert Main loop with [processx](https://github.com/r-lib/processx)

