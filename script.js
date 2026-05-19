// ===========================
// Reading Progress Bar
// ===========================
const progressBar = document.getElementById('reading-progress');
if (progressBar) {
  window.addEventListener('scroll', () => {
    const total = document.documentElement.scrollHeight - window.innerHeight;
    progressBar.style.width = (window.scrollY / total * 100) + '%';
  });
}

// ===========================
// Star Canvas Background
// ===========================
(function initStars() {
  const canvas = document.getElementById('star-canvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  let stars = [], W, H;

  function resize() {
    W = canvas.width  = window.innerWidth;
    H = canvas.height = window.innerHeight;
  }

  function createStars() {
    stars = [];
    const count = Math.floor((W * H) / 8000);
    for (let i = 0; i < count; i++) {
      stars.push({
        x: Math.random() * W,
        y: Math.random() * H,
        r: Math.random() * 1.4 + 0.2,
        alpha: Math.random(),
        speed: Math.random() * 0.3 + 0.05,
        twinkle: Math.random() * 0.02 + 0.005,
        dir: Math.random() > 0.5 ? 1 : -1
      });
    }
  }

  function draw() {
    ctx.clearRect(0, 0, W, H);
    stars.forEach(s => {
      s.alpha += s.twinkle * s.dir;
      if (s.alpha >= 1 || s.alpha <= 0.05) s.dir *= -1;
      ctx.beginPath();
      ctx.arc(s.x, s.y, s.r, 0, Math.PI * 2);
      ctx.fillStyle = `rgba(200, 220, 255, ${s.alpha})`;
      ctx.fill();
    });
    requestAnimationFrame(draw);
  }

  resize();
  createStars();
  draw();
  window.addEventListener('resize', () => { resize(); createStars(); });
})();

// ===========================
// Scroll Reveal (IntersectionObserver)
// ===========================
(function initReveal() {
  const selectors = '.reveal, .reveal-left, .reveal-right, .reveal-scale';
  const els = document.querySelectorAll(selectors);
  if (!els.length) return;

  const obs = new IntersectionObserver((entries) => {
    entries.forEach(e => {
      if (e.isIntersecting) {
        e.target.classList.add('visible');
        obs.unobserve(e.target);
      }
    });
  }, { threshold: 0.12 });

  els.forEach(el => obs.observe(el));
})();

// ===========================
// Animated Counters
// ===========================
function animateCounter(el) {
  const target = parseFloat(el.dataset.target);
  const suffix = el.dataset.suffix || '';
  const duration = 1800;
  const start = performance.now();

  function update(now) {
    const t = Math.min((now - start) / duration, 1);
    const ease = 1 - Math.pow(1 - t, 3);
    const val = target % 1 === 0
      ? Math.round(ease * target)
      : (ease * target).toFixed(1);
    el.textContent = val + suffix;
    if (t < 1) requestAnimationFrame(update);
  }
  requestAnimationFrame(update);
}

(function initCounters() {
  const counters = document.querySelectorAll('.stat-number[data-target]');
  if (!counters.length) return;

  const obs = new IntersectionObserver(entries => {
    entries.forEach(e => {
      if (e.isIntersecting) {
        animateCounter(e.target);
        obs.unobserve(e.target);
      }
    });
  }, { threshold: 0.5 });

  counters.forEach(c => obs.observe(c));
})();

// ===========================
// Typewriter Effect
// ===========================
function typeWriter(el, text, speed = 55) {
  el.textContent = '';
  el.classList.add('typewriter-cursor');
  let i = 0;
  const timer = setInterval(() => {
    el.textContent += text[i++];
    if (i >= text.length) {
      clearInterval(timer);
      el.classList.remove('typewriter-cursor');
    }
  }, speed);
}

const twEl = document.getElementById('typewriter-title');
if (twEl) {
  const text = twEl.dataset.text || twEl.textContent;
  twEl.textContent = '';
  setTimeout(() => typeWriter(twEl, text, 50), 400);
}

// ===========================
// Active Sidebar Link on Scroll
// ===========================
(function initSidebarActive() {
  const links = document.querySelectorAll('.sidebar a[href^="#"]');
  if (!links.length) return;

  const sections = Array.from(links).map(l => document.querySelector(l.getAttribute('href'))).filter(Boolean);

  const obs = new IntersectionObserver(entries => {
    entries.forEach(e => {
      if (e.isIntersecting) {
        links.forEach(l => l.classList.remove('active'));
        const link = document.querySelector(`.sidebar a[href="#${e.target.id}"]`);
        if (link) link.classList.add('active');
      }
    });
  }, { rootMargin: '-30% 0px -60% 0px' });

  sections.forEach(s => obs.observe(s));
})();

// ===========================
// Year in footer
// ===========================
const yearEl = document.getElementById('year');
if (yearEl) yearEl.textContent = new Date().getFullYear();

// ===========================
// Utility helpers
// ===========================
const R_UNIVERSAL = 8.31446261815324e3;
const fmt = (x, digits = 3) =>
  Number.isFinite(x) ? Number(x).toLocaleString(undefined, { maximumFractionDigits: digits }) : '—';
function getFloat(sel) {
  const el = document.querySelector(sel);
  if (!el) return NaN;
  const v = parseFloat(el.value);
  return Number.isFinite(v) ? v : NaN;
}
function MPaToPa(x) { return Number.isFinite(x) ? x * 1e6 : NaN; }

// ===========================
// Area Ratio helpers
// ===========================
function areaRatioFromMach(M, gamma) {
  const g = gamma;
  const term = (2 / (g + 1)) * (1 + (g - 1) / 2 * M * M);
  const pow = (g + 1) / (2 * (g - 1));
  return (1 / M) * Math.pow(term, pow);
}

function machFromAreaRatio(targetAR, gamma) {
  if (!Number.isFinite(targetAR) || !Number.isFinite(gamma)) return NaN;
  if (targetAR <= 1) {
    let lo = 1e-8, hi = 0.999999;
    for (let i = 0; i < 80; ++i) {
      const mid = 0.5 * (lo + hi);
      const f = areaRatioFromMach(mid, gamma) - targetAR;
      if (Math.abs(f) < 1e-9) return mid;
      if (areaRatioFromMach(lo, gamma) - targetAR > 0) { if (f > 0) lo = mid; else hi = mid; }
      else { if (f > 0) hi = mid; else lo = mid; }
    }
    return 0.5 * (lo + hi);
  } else {
    let lo = 1.0000001, hi = 200;
    const f_lo = areaRatioFromMach(lo, gamma) - targetAR;
    let fhi = areaRatioFromMach(hi, gamma) - targetAR;
    let attempts = 0;
    while (f_lo * fhi > 0 && attempts++ < 10) { hi *= 2; fhi = areaRatioFromMach(hi, gamma) - targetAR; }
    if (f_lo * fhi > 0) return NaN;
    for (let i = 0; i < 80; ++i) {
      const mid = 0.5 * (lo + hi);
      const fmid = areaRatioFromMach(mid, gamma) - targetAR;
      if (Math.abs(fmid) < 1e-10) return mid;
      if ((areaRatioFromMach(lo, gamma) - targetAR) * fmid <= 0) hi = mid; else lo = mid;
    }
    return 0.5 * (lo + hi);
  }
}

// ===========================
// Isentropic Rocket Performance
// ===========================
function calculateIsPerf() {
  const expansionRatio = getFloat('#expansionRatio');
  const gamma = getFloat('#specificHeatRatio');
  const Pc_MPa = getFloat('#chamberPressure');
  const T0 = getFloat('#chamberTemperature');
  const ambientPressure_MPa = getFloat('#ambientPressure');
  const molarMass = getFloat('#molarMass');
  const At = getFloat('#throatArea');
  const P0 = Number.isFinite(Pc_MPa) ? MPaToPa(Pc_MPa) : NaN;
  const Pa = Number.isFinite(ambientPressure_MPa) ? MPaToPa(ambientPressure_MPa) : NaN;

  document.querySelectorAll('.results td').forEach(cell => cell.textContent = '—');

  if (!Number.isFinite(expansionRatio) || !Number.isFinite(gamma)) {
    alert('Please provide at least the expansion ratio and specific heat ratio.');
    return;
  }

  const Me = machFromAreaRatio(expansionRatio, gamma);
  const Te_div_T0 = 1 / (1 + (gamma - 1) / 2 * Me * Me);
  const Pe_div_P0 = Math.pow(1 + (gamma - 1) / 2 * Me * Me, -gamma / (gamma - 1));
  const Te = Number.isFinite(T0) ? T0 * Te_div_T0 : NaN;
  const Pe = Number.isFinite(P0) ? P0 * Pe_div_P0 : NaN;

  let Rspec = NaN;
  if (Number.isFinite(molarMass)) Rspec = R_UNIVERSAL / molarMass;
  const Ve = Number.isFinite(Rspec) && Number.isFinite(Te) ? Me * Math.sqrt(gamma * Rspec * Te) : NaN;

  const Pe_div_P0_term = Math.pow(Pe_div_P0, (gamma - 1) / gamma);
  let Cf = NaN;
  if (Number.isFinite(P0)) {
    const Pa_for_Cf = Number.isFinite(Pa) ? Pa : 0;
    const term1 = Math.sqrt((2 * Math.pow(gamma, 2)) / (gamma - 1) * Math.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1)) * (1 - Pe_div_P0_term));
    const term2 = ((Pe - Pa_for_Cf) / P0) * expansionRatio;
    Cf = term1 + term2;
  }

  let mdot = NaN;
  if (Number.isFinite(At) && Number.isFinite(Rspec) && Number.isFinite(T0) && Number.isFinite(P0)) {
    const chokeFactor = Math.pow(2 / (gamma + 1), (gamma + 1) / (2 * (gamma - 1)));
    mdot = P0 * At * Math.sqrt(gamma / (Rspec * T0)) * chokeFactor;
  }

  let thrust = NaN, thrustVac = NaN;
  if (Number.isFinite(Cf) && Number.isFinite(P0) && Number.isFinite(At)) {
    thrust = Cf * P0 * At;
    const t1v = Math.sqrt((2 * Math.pow(gamma, 2)) / (gamma - 1) * Math.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1)) * (1 - Pe_div_P0_term));
    const t2v = (Pe / P0) * expansionRatio;
    thrustVac = (t1v + t2v) * P0 * At;
  }

  const g0 = 9.80665;
  const Isp = Number.isFinite(thrust) && Number.isFinite(mdot) && mdot > 0 ? thrust / (mdot * g0) : NaN;

  document.getElementById('vExit').textContent   = fmt(Ve);
  document.getElementById('pExit').textContent   = fmt(Pe);
  document.getElementById('tExit').textContent   = fmt(Te);
  document.getElementById('machExit').textContent = fmt(Me, 4);
  document.getElementById('pRatio').textContent  = fmt(Pe_div_P0, 5);
  document.getElementById('tRatio').textContent  = fmt(Te_div_T0, 5);
  document.getElementById('cf').textContent      = fmt(Cf, 4);
  document.getElementById('thrust').textContent  = fmt(thrust);
  document.getElementById('thrustVac').textContent = fmt(thrustVac);
  document.getElementById('isp').textContent     = fmt(Isp, 3);
  document.getElementById('mdot').textContent    = fmt(mdot);
}
const computeBtn = document.getElementById('compute');
if (computeBtn) computeBtn.addEventListener('click', calculateIsPerf);

// ===========================
// Isentropic Flow Table
// ===========================
function ratiosFromMach(M, g) {
  const T0T = 1 + ((g - 1) / 2) * M * M;
  const p0p = Math.pow(T0T, g / (g - 1));
  const rho0rho = Math.pow(T0T, 1 / (g - 1));
  const AAstar = areaRatioFromMach(M, g);
  const TTstar = 1 / T0T * (1 + (g - 1) / 2);
  const ppstar = (1 / p0p) / Math.pow(2 / (g + 1), g / (g - 1));
  const rhorhostar = 1 / rho0rho / Math.pow(2 / (g + 1), 1 / (g - 1));
  return { T0T, p0p, rho0rho, AAstar, TTstar, ppstar, rhorhostar };
}
function machFromT0T(val, g) { return Math.sqrt(2 * (val - 1) / (g - 1)); }
function machFromP0P(val, g) { return machFromT0T(Math.pow(val, (g - 1) / g), g); }
function machFromRho0Rho(val, g) { return machFromT0T(Math.pow(val, g - 1), g); }
function machFromARHelper(targetAR, g) {
  let lo = 1e-6, hi = 20;
  for (let i = 0; i < 100; i++) {
    const mid = 0.5 * (lo + hi);
    const f = areaRatioFromMach(mid, g) - targetAR;
    if (Math.abs(f) < 1e-6) return mid;
    if (f > 0) hi = mid; else lo = mid;
  }
  return 0.5 * (lo + hi);
}

function fillIsentropicTable() {
  const g = parseFloat(document.getElementById('gamma').value);
  const param = document.getElementById('knownParam').value;
  const val = parseFloat(document.getElementById('paramValue').value);
  if (!Number.isFinite(g) || !Number.isFinite(val)) { alert('Please enter valid inputs.'); return; }
  let M;
  if (param === 'Ma')       M = val;
  if (param === 'T0T')      M = machFromT0T(val, g);
  if (param === 'p0p')      M = machFromP0P(val, g);
  if (param === 'rho0rho')  M = machFromRho0Rho(val, g);
  if (param === 'AAstar')   M = machFromARHelper(val, g);
  if (!Number.isFinite(M)) { alert('Could not determine Mach number.'); return; }
  const r = ratiosFromMach(M, g);
  document.getElementById('res-gamma').textContent    = fmt(g);
  document.getElementById('res-Ma').textContent       = fmt(M);
  document.getElementById('res-T0T').textContent      = fmt(r.T0T);
  document.getElementById('res-p0p').textContent      = fmt(r.p0p);
  document.getElementById('res-rho0rho').textContent  = fmt(r.rho0rho);
  document.getElementById('res-AAstar').textContent   = fmt(r.AAstar);
  document.getElementById('res-TTstar').textContent   = fmt(r.TTstar);
  document.getElementById('res-ppstar').textContent   = fmt(r.ppstar);
  document.getElementById('res-rhorhostar').textContent = fmt(r.rhorhostar);
}
const computeISBtn = document.getElementById('computeIS');
if (computeISBtn) computeISBtn.addEventListener('click', fillIsentropicTable);

// ===========================
// Drag / Lift Calculator
// ===========================
function computeForces() {
  const speed = parseFloat(document.getElementById('speed').value);
  const rho   = parseFloat(document.getElementById('rho').value);
  const area  = parseFloat(document.getElementById('area').value);
  const cd    = parseFloat(document.getElementById('cd').value);
  const cl    = parseFloat(document.getElementById('cl').value);
  if (isNaN(speed) || isNaN(rho) || isNaN(area)) {
    ['q','drag','lift'].forEach(id => document.getElementById(id).textContent = '—');
    return;
  }
  const q = 0.5 * rho * speed * speed;
  document.getElementById('q').textContent    = q.toFixed(2);
  document.getElementById('drag').textContent = isNaN(cd) ? '—' : (q * area * cd).toFixed(2);
  document.getElementById('lift').textContent = isNaN(cl) ? '—' : (q * area * cl).toFixed(2);
}
const dragBtn = document.getElementById('computedrag');
if (dragBtn) dragBtn.addEventListener('click', computeForces);

// ===========================
// Oblique Shock Calculator
// ===========================
function computeObliqueShock() {
  const M1       = parseFloat(document.getElementById('mach1').value);
  const gamma    = parseFloat(document.getElementById('gamma_os').value);
  const betaDeg  = parseFloat(document.getElementById('beta').value);
  const ids = ['mach2','theta','pRatio_os','rhoRatio_os','tRatio_os','p0Ratio_os'];
  if (isNaN(M1) || isNaN(gamma) || isNaN(betaDeg) || M1 <= 1) {
    ids.forEach(id => document.getElementById(id).textContent = '—'); return;
  }
  const betaRad = betaDeg * (Math.PI / 180);
  const Mn1 = M1 * Math.sin(betaRad);
  if (Mn1 <= 1) {
    ids.forEach(id => document.getElementById(id).textContent = '—');
    alert('Normal component of Mach must be > 1 for a shock.'); return;
  }
  const pRatio   = ((2 * gamma) / (gamma + 1)) * Mn1 * Mn1 - ((gamma - 1) / (gamma + 1));
  const rhoRatio = ((gamma + 1) * Mn1 * Mn1) / ((gamma - 1) * Mn1 * Mn1 + 2);
  const tRatio   = pRatio / rhoRatio;
  const Mn2sq    = (Mn1 * Mn1 * (gamma - 1) + 2) / (2 * gamma * Mn1 * Mn1 - (gamma - 1));
  const M2       = Math.sqrt(Mn2sq / Math.sin(betaRad) ** 2 + (1 / Math.tan(betaRad)) ** 2);
  const thetaRad = Math.atan(2 * (M1 * M1 * Math.sin(betaRad) ** 2 - 1) * Math.cos(betaRad) / (M1 * M1 * (gamma + Math.cos(2 * betaRad)) + 2));
  const t1 = Math.pow(((gamma + 1) * (M1 * Math.sin(betaRad)) ** 2) / ((gamma - 1) * (M1 * Math.sin(betaRad)) ** 2 + 2), gamma / (gamma - 1));
  const t2 = Math.pow((gamma + 1) / (2 * gamma * (M1 * Math.sin(betaRad)) ** 2 - (gamma - 1)), 1 / (gamma - 1));
  document.getElementById('mach2').textContent    = M2.toFixed(3);
  document.getElementById('theta').textContent    = (thetaRad * 180 / Math.PI).toFixed(2);
  document.getElementById('pRatio_os').textContent  = pRatio.toFixed(3);
  document.getElementById('rhoRatio_os').textContent = rhoRatio.toFixed(3);
  document.getElementById('tRatio_os').textContent  = tRatio.toFixed(3);
  document.getElementById('p0Ratio_os').textContent = (t1 * t2).toFixed(3);
}
const obliqueBtn = document.getElementById('computeOblique');
if (obliqueBtn) obliqueBtn.addEventListener('click', computeObliqueShock);

(function initRocket() {
if (!document.querySelector('.orbit-deco') || window.innerWidth <= 768) return;  const c = document.createElement('canvas');
  c.style.cssText = 'position:fixed;top:0;left:0;width:100%;height:100%;pointer-events:none;z-index:2;';
  document.body.appendChild(c);
  const ctx = c.getContext('2d');
  let W, H, particles = [];

  function resize() { W = c.width = window.innerWidth; H = c.height = window.innerHeight; }
  resize();
  window.addEventListener('resize', resize);

  function draw() {
    ctx.clearRect(0, 0, W, H);
    const startY      = 0.535;
    const scrollSpeed = 0.12;
    const maxScroll   = document.documentElement.scrollHeight - window.innerHeight;
    const progress    = Math.min(window.scrollY / (maxScroll * scrollSpeed), 1);
    const y = H * startY - progress * (H * startY + 150);
    const x = W / 2;

    // Flame particles (adjusted 'y' slightly to fit the new engine nozzle nozzle)
    if (window.scrollY > 10 && progress < 1 && Math.random() > 0.3)
      particles.push({ x: x + (Math.random()-0.5)*12, y: y+52, vx: (Math.random()-0.5)*1.5, vy: Math.random()*3+2, life: 1 });
    
    particles = particles.filter(p => p.life > 0);
    particles.forEach(p => {
      p.x += p.vx; p.y += p.vy; p.life -= 0.04;
      const g = ctx.createRadialGradient(p.x,p.y,0,p.x,p.y,12);
      g.addColorStop(0,`rgba(0,245,212,${p.life})`);
      g.addColorStop(1,`rgba(0, 13, 255, 0)`);
      ctx.fillStyle = g;
      ctx.beginPath(); ctx.arc(p.x,p.y,12,0,Math.PI*2); ctx.fill();
    });

    if (y > -150) {
      // 1. Create lighting gradients for a 3D cylindrical effect
      const bodyGrad = ctx.createLinearGradient(x - 20, 0, x + 20, 0);
      bodyGrad.addColorStop(0, '#ffffff');
      bodyGrad.addColorStop(0.3, '#e7eaf6');
      bodyGrad.addColorStop(1, '#4e4f53');

      const accentGrad = ctx.createLinearGradient(x - 20, 0, x + 20, 0);
      accentGrad.addColorStop(0, '#19ae9a');
      accentGrad.addColorStop(1, '#025f74');

      // 2. Engine Nozzle (Dark metallic trapezoid at the base)
      ctx.fillStyle = '#3d405b';
      ctx.beginPath();
      ctx.moveTo(x - 10, y + 40);
      ctx.lineTo(x + 10, y + 40);
      ctx.lineTo(x + 14, y + 50);
      ctx.lineTo(x - 14, y + 50);
      ctx.closePath();
      ctx.fill();

      // 3. Wings / Fins (Swept-back aerodynamic curves)
      ctx.fillStyle = accentGrad;
      // Left Fin
      ctx.beginPath();
      ctx.moveTo(x - 16, y + 10);
      ctx.bezierCurveTo(x - 35, y + 15, x - 42, y + 45, x - 38, y + 55);
      ctx.lineTo(x - 16, y + 42);
      ctx.closePath();
      ctx.fill();

      // Right Fin
      ctx.beginPath();
      ctx.moveTo(x + 16, y + 10);
      ctx.bezierCurveTo(x + 35, y + 15, x + 42, y + 45, x + 38, y + 55);
      ctx.lineTo(x + 16, y + 42);
      ctx.closePath();
      ctx.fill();

      // 4. Main Fuselage Body (Sleek capsule structure)
      ctx.fillStyle = bodyGrad;
      ctx.beginPath();
      ctx.moveTo(x - 16, y - 20);
      ctx.lineTo(x + 16, y - 20);
      ctx.quadraticCurveTo(x + 19, y + 10, x + 16, y + 40);
      ctx.lineTo(x - 16, y + 40);
      ctx.quadraticCurveTo(x - 19, y + 10, x - 16, y - 20);
      ctx.closePath();
      ctx.fill();

      // 5. Nose Cone (Smoothly tapered point)
      ctx.fillStyle = accentGrad;
      ctx.beginPath();
      ctx.moveTo(x - 16, y - 20);
      ctx.quadraticCurveTo(x - 14, y - 52, x, y - 68); // Soft curve to the tip
      ctx.quadraticCurveTo(x + 14, y - 52, x + 16, y - 20);
      ctx.closePath();
      ctx.fill();

      // 6. Window / Porthole (Circular with a metallic rim and a gloss layer)
      ctx.fillStyle = '#b5bad0'; // Outer Rim
      ctx.beginPath(); ctx.arc(x, y + 5, 11, 0, Math.PI * 2); ctx.fill();

      ctx.fillStyle = '#4db8ff'; // Glass
      ctx.beginPath(); ctx.arc(x, y + 5, 8, 0, Math.PI * 2); ctx.fill();
      
      ctx.fillStyle = 'rgba(255, 255, 255, 0.4)'; // Glass Reflection Shine
      ctx.beginPath();
      ctx.arc(x, y + 5, 8, Math.PI * 1.2, Math.PI * 1.7);
      ctx.lineTo(x, y + 5);
      ctx.closePath();
      ctx.fill();
    }
    requestAnimationFrame(draw);
  }
  draw();
})();