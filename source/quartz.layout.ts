import { PageLayout, SharedLayout } from "./quartz/cfg"
import * as Component from "./quartz/components"
import { UtterancesScript } from "./quartz/components/UtterancesScript"

// components shared across all pages
export const sharedPageComponents: SharedLayout = {
  head: Component.Head({ extraScripts: [UtterancesScript()] }),
  header: [],
  afterBody: [
    Component.Comments({
      provider: 'giscus',
      options: {
        repo: 'LabOnoM/DK.BeesGO',
        repoId: 'R_kgDOL6VQ_Q',
        category: 'Announcements',
        categoryId: 'DIC_kwDOL6VQ_c4CfaND',
        mapping: 'pathname',
        inputPosition: 'top',
        reactionsEnabled: true,
        strict: true,
        themeUrl: 'https://giscus.app/themes/', 
        lightTheme: 'noborder_light', 
        darkTheme: 'noborder_dark',
      }
    }),
  ],
  footer: Component.Footer({
    links: {
      "BSGOU": "https://www.bs-gou.com/",
      GitHub: "https://github.com/jackyzha0/quartz",
      "Discord Community": "https://discord.gg/cRFFHYye7t",
    },
  }),
}

// components for pages that display a single page (e.g. a single note)
export const defaultContentPageLayout: PageLayout = {
  beforeBody: [
    Component.Breadcrumbs(),
    Component.ArticleTitle(),
    Component.ContentMeta(),
    Component.TagList(),
  ],
  left: [
    Component.PageTitle(),
    Component.MobileOnly(Component.Spacer()),
    Component.Search(),
    Component.Darkmode(),
    Component.DesktopOnly(Component.Explorer()),
  ],
  right: [
    Component.Graph(),
    Component.DesktopOnly(Component.TableOfContents()),
    Component.Backlinks(),
  ],
}

// components for pages that display lists of pages  (e.g. tags or folders)
export const defaultListPageLayout: PageLayout = {
  beforeBody: [Component.Breadcrumbs(), Component.ArticleTitle(), Component.ContentMeta()],
  left: [
    Component.PageTitle(),
    Component.MobileOnly(Component.Spacer()),
    Component.Search(),
    Component.Darkmode(),
    Component.DesktopOnly(Component.Explorer()),
  ],
  right: [],
}
