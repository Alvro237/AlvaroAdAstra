/* ===========================================================
   Idioma de la web (ES / EN)
   - Cada texto existe dos veces: data-lang="en" y data-lang="es".
     style.css oculta el que no toca según <html lang="...">.
   - Atributos traducibles: data-es-placeholder, data-es-alt,
     data-es-aria-label, data-es-content, data-es-title...
     (el valor en inglés es el propio atributo). data-es-text cambia
     el texto del elemento (para <option>, que no admite <span>).
   - <title data-es="..."> traduce el título de la pestaña.
   - Prioridad: ?lang=es|en en la URL > elección guardada > idioma del navegador.
   - Al cambiar dispara el evento "sitelangchange" (lo usa el mapa).
   Se carga en <head> SIN defer para que no parpadee el idioma equivocado.
   =========================================================== */
(function () {
  'use strict';
  var KEY = 'site-lang';
  var html = document.documentElement;

  function valida(l) { return l === 'es' || l === 'en' ? l : null; }
  function leer() {
    try { return valida(localStorage.getItem(KEY)); } catch (e) { return null; }
  }
  function guardar(l) {
    try { localStorage.setItem(KEY, l); } catch (e) { /* modo privado: no pasa nada */ }
  }

  var deUrl = null;
  try { deUrl = valida(new URLSearchParams(location.search).get('lang')); } catch (e) { /* navegador antiguo */ }
  var nav = (navigator.languages && navigator.languages[0]) || navigator.language || '';
  var lang = deUrl || leer() || (/^es\b/i.test(nav) ? 'es' : 'en');
  if (deUrl) guardar(deUrl);
  html.lang = lang;

  function traducirAtributos(raiz) {
    var els = (raiz || document).querySelectorAll('*');
    for (var i = 0; i < els.length; i++) {
      var el = els[i], attrs = el.attributes;
      for (var j = 0; j < attrs.length; j++) {
        var n = attrs[j].name;
        if (n.indexOf('data-es-') !== 0) continue;
        var destino = n.slice(8);
        var copiaEn = 'data-en-' + destino;
        if (destino === 'text') {
          // Para elementos que no admiten <span> dentro (p. ej. <option>): se cambia el texto.
          if (!el.hasAttribute(copiaEn)) el.setAttribute(copiaEn, el.textContent);
          el.textContent = lang === 'es' ? attrs[j].value : el.getAttribute(copiaEn);
          continue;
        }
        if (!el.hasAttribute(copiaEn)) el.setAttribute(copiaEn, el.getAttribute(destino) || '');
        el.setAttribute(destino, lang === 'es' ? attrs[j].value : el.getAttribute(copiaEn));
      }
    }
    var t = document.querySelector('title[data-es]');
    if (t) {
      if (!t.hasAttribute('data-en')) t.setAttribute('data-en', t.textContent);
      document.title = lang === 'es' ? t.getAttribute('data-es') : t.getAttribute('data-en');
    }
  }

  function pintarBoton() {
    var b = document.getElementById('lang-switch');
    if (!b) return;
    b.setAttribute('aria-label', lang === 'es' ? 'Switch to English' : 'Cambiar a español');
    b.title = b.getAttribute('aria-label');
    b.querySelector('[data-l="es"]').classList.toggle('on', lang === 'es');
    b.querySelector('[data-l="en"]').classList.toggle('on', lang === 'en');
  }

  function poner(l) {
    l = valida(l);
    if (!l || l === lang) return;
    lang = l;
    html.lang = l;
    guardar(l);
    traducirAtributos();
    pintarBoton();
    var ev;
    try { ev = new CustomEvent('sitelangchange', { detail: { lang: l } }); }
    catch (e) { ev = document.createEvent('CustomEvent'); ev.initCustomEvent('sitelangchange', false, false, { lang: l }); }
    document.dispatchEvent(ev);
  }

  window.SiteLang = {
    get: function () { return lang; },
    set: poner,
    /** Devuelve el texto en el idioma activo: SiteLang.t({ en: '...', es: '...' }) */
    t: function (o) { return o && (o[lang] != null ? o[lang] : o.en); }
  };

  function iniciar() {
    var navEl = document.querySelector('header nav');
    if (navEl && !document.getElementById('lang-switch')) {
      var b = document.createElement('button');
      b.type = 'button';
      b.id = 'lang-switch';
      b.className = 'lang-switch';
      b.innerHTML = '<span data-l="es">ES</span><span data-l="en">EN</span>';
      b.addEventListener('click', function () { poner(lang === 'es' ? 'en' : 'es'); });
      navEl.appendChild(b);
    }
    traducirAtributos();
    pintarBoton();
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', iniciar);
  else iniciar();
})();
