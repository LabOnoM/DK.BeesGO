document.addEventListener("nav", () => {
  const openBtn = document.getElementById("open-subscribe-popup") as HTMLButtonElement | null
  const popup = document.getElementById("subscribe-popup") as HTMLElement | null
  const closeBtn = document.getElementById("close-subscribe-popup") as HTMLButtonElement | null
  if (openBtn && popup && closeBtn) {
    const open = () => {
      popup.style.display = "block"
    }
    const close = () => {
      popup.style.display = "none"
    }
    const outside = (e: MouseEvent) => {
      if (e.target === popup) {
        popup.style.display = "none"
      }
    }
    openBtn.addEventListener("click", open)
    closeBtn.addEventListener("click", close)
    popup.addEventListener("click", outside)
    window.addCleanup(() => {
      openBtn.removeEventListener("click", open)
      closeBtn.removeEventListener("click", close)
      popup.removeEventListener("click", outside)
    })
  }
})
