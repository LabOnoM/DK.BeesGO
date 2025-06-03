window.addEventListener("DOMContentLoaded", function() {
  // Only add comments section to content pages, not lists/tags
  // Check for unique selector for your notes
  if (document.querySelector("main")) {
    const section = document.createElement('section');
    section.id = 'comments';
    section.style.marginTop = '3em';
    document.querySelector("main").appendChild(section);

    const s = document.createElement('script');
    s.src = "https://utteranc.es/client.js";
    s.setAttribute("repo", "LabOnoM/DK.BeesGO");
    s.setAttribute("issue-term", "pathname");
    s.setAttribute("theme", "github-light");
    s.crossOrigin = "anonymous";
    s.async = true;
    section.appendChild(s);
  }
});
