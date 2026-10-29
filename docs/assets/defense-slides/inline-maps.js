/* Lazy in-slide atlas frames, static fallbacks and an in-page enlargement. */
(function () {
  "use strict";
  const es = /^es\b/i.test(document.documentElement.lang || "en");
  const tr = (en, spanish) => es ? spanish : en;
  const entries = Array.from(document.querySelectorAll(".inline-atlas")).map(root => ({
    root, frame: root.querySelector("iframe"),
    status: root.querySelector(".map-load-status"),
    toggle: root.querySelector(".map-static-toggle"),
    gene: root.dataset.defaultGene, preferStatic: false, loaded: false, expanded: false
  }));
  if (!entries.length) return;
  document.querySelectorAll("a[data-map-slide]").forEach(link => {
    link.addEventListener("click", event => {
      event.preventDefault();
      event.stopPropagation();
      const target = document.getElementById(link.dataset.mapSlide);
      const index = Array.from(document.querySelectorAll(".reveal .slides > section")).indexOf(target);
      if (index >= 0 && window.Reveal && window.Reveal.slide) window.Reveal.slide(index);
      else if (target) window.location.hash = "/" + target.id;
    });
  });

  const dialog = document.createElement("dialog");
  dialog.className = "inline-map-dialog";
  dialog.setAttribute("aria-label", tr("Enlarged interactive brain map", "Mapa cerebral interactivo ampliado"));
  const bar = document.createElement("div");
  bar.className = "inline-map-dialog-bar";
  const heading = document.createElement("h2");
  const close = document.createElement("button");
  close.type = "button";
  close.textContent = tr("Return to slide", "Volver a la diapositiva");
  close.addEventListener("click", () => dialog.close());
  bar.append(heading, close);
  const largeFrame = document.createElement("iframe");
  largeFrame.className = "inline-map-dialog-frame";
  largeFrame.title = tr("Interactive reference and PM2.5 scenario", "Referencia interactiva y escenario PM2.5");
  dialog.append(bar, largeFrame);
  document.body.appendChild(dialog);
  let enlarged = null;
  let instance = 0;

  function activeSlide() {
    return window.Reveal && window.Reveal.getCurrentSlide
      ? window.Reveal.getCurrentSlide() : document.querySelector(".reveal .slides > section.present");
  }
  function isActive(entry) { return entry.root.closest("section") === activeSlide(); }
  function mapURL(entry) {
    const url = new URL(entry.frame.dataset.mapSrc, window.location.href);
    if (entry.gene) url.searchParams.set("defaultgene", entry.gene);
    url.searchParams.set("instance", entry.instance || "");
    return url;
  }
  function release(entry) {
    if (entry.loaded) entry.frame.src = "about:blank";
    entry.loaded = false;
    entry.root.classList.remove("map-live");
    clearTimeout(entry.timer);
  }
  function resize(frame, origin) {
    if (frame.contentWindow) frame.contentWindow.postMessage({type:"omics-slide-map-resize"}, origin === "null" ? "*" : origin);
  }
  function load(entry) {
    if (entry.loaded || entry.preferStatic || entry.expanded || !isActive(entry)) return;
    const url = mapURL(entry);
    entry.instance = String(++instance);
    url.searchParams.set("instance", entry.instance);
    entry.loaded = true;
    entry.status.textContent = tr("Loading interactive map…", "Cargando mapa interactivo…");
    entry.frame.src = url.href;
    entry.timer = setTimeout(() => {
      if (entry.loaded && !entry.root.classList.contains("map-live")) {
        release(entry);
        entry.preferStatic = true;
        entry.toggle.textContent = tr("Activate 3D", "Activar 3D");
        entry.toggle.setAttribute("aria-pressed", "true");
        entry.status.textContent = tr("The map has not loaded. Reactivate 3D to retry.", "El mapa no ha cargado. Reactive 3D para reintentar.");
      }
    }, 30000);
  }
  function focusSlides(entry) {
    const slide = entry.root.closest("section");
    slide.tabIndex = -1;
    slide.focus({preventScroll:true});
  }
  function sync() {
    if (enlarged && !isActive(enlarged)) dialog.close();
    entries.forEach(entry => {
      if (isActive(entry)) load(entry);
      else release(entry);
    });
  }
  entries.forEach(entry => {
    entry.toggle.addEventListener("click", () => {
      entry.preferStatic = !entry.preferStatic;
      entry.toggle.textContent = entry.preferStatic ? tr("Activate 3D", "Activar 3D") : tr("Static image", "Imagen estática");
      entry.toggle.setAttribute("aria-pressed", String(entry.preferStatic));
      if (entry.preferStatic) {
        release(entry);
        entry.status.textContent = tr("Static example shown", "Ejemplo estático visible");
      } else load(entry);
    });
    entry.root.querySelector(".map-focus-slides").addEventListener("click", () => focusSlides(entry));
    entry.root.querySelector(".map-enlarge").addEventListener("click", () => {
      if (enlarged) dialog.close();
      enlarged = entry;
      entry.expanded = true;
      release(entry); // One atlas frame at a time avoids accumulating WebGL contexts.
      entry.instance = String(++instance);
      heading.textContent = tr("Interactive ", "Mapa interactivo · ") + (entry.root.dataset.region === "HY" ? tr("hypothalamus", "hipotálamo") : tr("cortex", "corteza"));
      dialog.showModal();
      largeFrame.src = mapURL(entry).href;
    });
  });
  dialog.addEventListener("keydown", event => event.stopPropagation());
  dialog.addEventListener("close", () => {
    largeFrame.src = "about:blank";
    const previous = enlarged;
    enlarged = null;
    if (previous) {
      previous.expanded = false;
      load(previous);
      focusSlides(previous);
    }
  });
  window.addEventListener("message", event => {
    const data = event.data;
    if (!data || data.type !== "omics-slide-map-state") return;
    const entry = entries.find(item => event.source === item.frame.contentWindow) ||
      (enlarged && event.source === largeFrame.contentWindow ? enlarged : null);
    if (!entry) return;
    const expected = mapURL(entry);
    const validOrigin = expected.protocol === "file:"
      ? window.location.protocol === "file:" && event.origin === "null"
      : event.origin === expected.origin;
    if (!validOrigin || data.region !== entry.root.dataset.region || data.instance !== entry.instance) return;
    const selected = data.status === "loading" ? data.requestedGene : data.gene;
    if (typeof selected === "string" && /^[A-Za-z0-9_.-]{1,80}$/.test(selected)) entry.gene = selected;
    if (data.status === "ready") {
      clearTimeout(entry.timer);
      if (event.source === entry.frame.contentWindow && entry.loaded && isActive(entry) && !entry.expanded && !entry.preferStatic) entry.root.classList.add("map-live");
      entry.status.textContent = tr("Interactive map · ", "Mapa interactivo · ") + entry.gene;
      requestAnimationFrame(() => resize(event.source === largeFrame.contentWindow ? largeFrame : entry.frame, expected.origin));
    } else if (data.status === "error") {
      entry.status.textContent = tr("Map unavailable. Use the static image or reactivate 3D.", "Mapa no disponible. Use la imagen estática o reactive 3D.");
      if (event.source === entry.frame.contentWindow) {
        release(entry);
        entry.preferStatic = true;
        entry.toggle.textContent = tr("Activate 3D", "Activar 3D");
        entry.toggle.setAttribute("aria-pressed", "true");
      }
    }
  });
  window.addEventListener("hashchange", sync);
  window.addEventListener("resize", () => entries.forEach(entry => {
    if (entry.loaded) resize(entry.frame, mapURL(entry).origin);
  }));
  if (window.Reveal && window.Reveal.on) {
    window.Reveal.on("ready", sync);
    window.Reveal.on("slidechanged", sync);
  }
  sync();
}());
