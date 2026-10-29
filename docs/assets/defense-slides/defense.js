/* Language navigation and an offline-readable notes drawer. */
(function () {
  "use strict";
  var language = (document.documentElement.lang || "en").toLowerCase();
  var spanish = language.indexOf("es") === 0;
  var current = new URL(window.location.href);
  var alternate = document.querySelector('link[rel="alternate"][hreflang="' + (spanish ? "en" : "es") + '"]');
  var destination;
  if (alternate && alternate.href) {
    destination = new URL(alternate.href, current);
  } else {
    destination = new URL(current);
    if (spanish) {
      destination.pathname = destination.pathname.replace(/\/es\/defense\.html$/, "/defense.html");
    } else {
      destination.pathname = destination.pathname.replace(/\/defense\.html$/, "/es/defense.html");
    }
  }
  // No speculative language URL is exposed for another filename or route.
  if (destination.pathname === current.pathname) return;
  destination.search = current.search;
  destination.hash = current.hash;
  var link = document.createElement("a");
  link.className = "defense-language-switch";
  link.href = destination.href;
  link.target = "_self";
  link.addEventListener("click", function (event) {
    event.preventDefault();
    window.location.assign(link.href);
  });
  link.hreflang = spanish ? "en" : "es";
  link.lang = spanish ? "en" : "es";
  link.textContent = spanish ? "English" : "Español";
  link.setAttribute("aria-label", spanish ? "Ver esta diapositiva en inglés" : "View this slide in Spanish");
  document.body.appendChild(link);
  var notesButton = document.createElement("button");
  notesButton.type = "button";
  notesButton.className = "defense-notes-toggle";
  notesButton.textContent = spanish ? "Notas" : "Notes";
  notesButton.setAttribute("aria-expanded", "false");
  notesButton.setAttribute("aria-controls", "defense-notes");
  var notesPanel = document.createElement("section");
  notesPanel.className = "defense-notes";
  notesPanel.id = "defense-notes";
  notesPanel.hidden = true;
  notesPanel.tabIndex = -1;
  notesPanel.setAttribute("role", "dialog");
  notesPanel.setAttribute("aria-modal", "false");
  notesPanel.setAttribute("aria-label", spanish ? "Notas de la diapositiva actual" : "Notes for the current slide");
  var closeButton = document.createElement("button");
  closeButton.type = "button";
  closeButton.className = "defense-notes-close";
  closeButton.textContent = spanish ? "Cerrar" : "Close";
  var notesHeading = document.createElement("h2");
  var notesBody = document.createElement("div");
  notesPanel.appendChild(closeButton);
  notesPanel.appendChild(notesHeading);
  notesPanel.appendChild(notesBody);
  document.body.appendChild(notesButton);
  document.body.appendChild(notesPanel);
  function activeSlide() {
    return window.Reveal && typeof window.Reveal.getCurrentSlide === "function"
      ? window.Reveal.getCurrentSlide()
      : document.querySelector(".reveal .slides section.present");
  }
  function updateNotes() {
    var active = activeSlide();
    var heading = active && active.querySelector("h1, h2, h3");
    var notes = active && active.querySelector("aside.notes");
    notesHeading.textContent = heading ? heading.textContent : (spanish ? "Notas" : "Notes");
    notesBody.innerHTML = notes ? notes.innerHTML : "";
    if (!notes) notesBody.textContent = spanish ? "Esta diapositiva no tiene notas adicionales." : "This slide has no additional notes.";
    notesPanel.scrollTop = 0;
  }
  function closeNotes() {
    notesPanel.hidden = true;
    notesButton.setAttribute("aria-expanded", "false");
    notesButton.focus();
  }
  notesButton.addEventListener("click", function () {
    if (!notesPanel.hidden) { closeNotes(); return; }
    updateNotes();
    notesPanel.hidden = false;
    notesButton.setAttribute("aria-expanded", "true");
    notesPanel.focus();
  });
  closeButton.addEventListener("click", closeNotes);
  notesPanel.addEventListener("keydown", function (event) {
    // Arrow/space keys scroll the readable notes rather than advancing the deck.
    event.stopPropagation();
    if (event.key === "Escape") { event.preventDefault(); closeNotes(); }
  });
  function sync() {
    destination.hash = window.location.hash;
    link.href = destination.href;
    var active = activeSlide();
    document.body.dataset.defenseDark = active && (active.classList.contains("section-slide") || active.id === "title-slide") ? "true" : "false";
    if (!notesPanel.hidden) updateNotes();
  }
  window.addEventListener("hashchange", sync);
  function bindReveal() {
    if (window.Reveal && typeof window.Reveal.on === "function") {
      window.Reveal.on("slidechanged", sync);
      window.Reveal.on("ready", sync);
    }
    sync();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", bindReveal, { once: true });
  else bindReveal();
}());
