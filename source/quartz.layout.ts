import { PageLayout, SharedLayout } from "./quartz/cfg"
import * as Component from "./quartz/components"

// components shared across all pages
export const sharedPageComponents: SharedLayout = {
  head: Component.Head(),
  header: [],
  afterBody: [
    Component.Comments({
      provider: 'giscus',
      options: {
        repo: 'LabOnoM/DK.BeesGO',
        repoId: 'R_kgDOO1GkSQ',
        category: 'General',
        categoryId: 'DIC_kwDOO1GkSc4Cq-H5',
        mapping: 'pathname',
        strict: false,
        inputPosition: 'top',
        reactionsEnabled: true,
        themeUrl: 'https://giscus.app/themes/', 
        lightTheme: 'noborder_light', 
        darkTheme: 'noborder_dark',
      }
    }),
  ],
  footer: Component.Footer({
    links: {
      "Home": "https://www.bs-gou.com/index.html",
      "Blog": "https://www.bs-gou.com/blog/index.html",
      "GitHub": "https://github.com/LabOnoM",
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
    Component.DesktopOnly(Component.Explorer({
      title: "Bees-GO!",
      filterFn: (node) => node.displayName.toLowerCase() !== "tags",
    })),
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
