// ---------- Utility helpers ----------
const R_UNIVERSAL = 8.31446261815324e3; // J / (kmol K)

const fmt = (x, digits = 3) =>
  Number.isFinite(x)
    ? Number(x).toLocaleString(undefined, { maximumFractionDigits: digits })
    : "—";

function getFloat(sel) {
  const el = document.querySelector(sel);
  if (!el) return NaN;
  const v = parseFloat(el.value);
  return Number.isFinite(v) ? v : NaN;
}

function MPaToPa(x) {
  return Number.isFinite(x) ? x * 1e6 : NaN;
}

// ---------- Area-Mach relation ----------
function areaRatioFromMach(M, g) {
  const term = (2 / (g + 1)) * (1 + ((g - 1) / 2) * M * M);
  const pow = (g + 1) / (2 * (g - 1));
  return (1 / M) * Math.pow(term, pow);
}

function machFromAreaRatio(targetAR, g) {
  if (!Number.isFinite(targetAR) || !Number.isFinite(g)) return NaN;
  if (targetAR <= 1) {
    // subsonic branch
    let lo = 1e-8,
      hi = 0.999999;
    for (let i = 0; i < 80; ++i) {
      const mid = 0.5 * (lo + hi);
      const f = areaRatioFromMach(mid, g) - targetAR;
      if (Math.abs(f) < 1e-9) return mid;
      if (areaRatioFromMach(lo, g) - targetAR > 0) {
        if (f > 0) lo = mid;
        else hi = mid;
      } else {
        if (f > 0) hi = mid;
        else lo = mid;
      }
    }
    return 0.5 * (lo + hi);
  } else {
    // supersonic branch
    let lo = 1.0000001,
      hi = 200;
    for (let i = 0; i < 80; ++i) {
      const mid = 0.5 * (lo + hi);
      const fmid = areaRatioFromMach(mid, g) - targetAR;
      if (Math.abs(fmid) < 1e-10) return mid;
      if ((areaRatioFromMach(lo, g) - targetAR) * fmid <= 0) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    return 0.5 * (lo + hi);
  }
}

// ---------- Main performance calculation ----------
function calculateIsPerf() {
  const expansionRatio = getFloat("#expansionRatio");
  const gamma = getFloat("#specificHeatRatio");
  let Pc_MPa = getFloat("#chamberPressure");
  const T0 = getFloat("#chamberTemperature");

  let ambientPressure_MPa = getFloat("#ambientPressure");
  const molarMass = getFloat("#molarMass");
  const At = getFloat("#throatArea");

  const P0 = Number.isFinite(Pc_MPa) ? MPaToPa(Pc_MPa) : NaN;
  const Pa = Number.isFinite(ambientPressure_MPa)
    ? MPaToPa(ambientPressure_MPa)
    : NaN;

  document
    .querySelectorAll(".results td")
    .forEach((cell) => (cell.textContent = "—"));

  if (!Number.isFinite(expansionRatio) || !Number.isFinite(gamma)) {
    alert("Please provide at least the expansion ratio and specific heat ratio.");
    return;
  }

  const Me = machFromAreaRatio(expansionRatio, gamma);

  const Te_div_T0 = 1 / (1 + ((gamma - 1) / 2) * Me * Me);
  const Pe_div_P0 = Math.pow(1 + ((gamma - 1) / 2) * Me * Me, -gamma / (gamma - 1));

  const Te = Number.isFinite(T0) ? T0 * Te_div_T0 : NaN;
  const Pe = Number.isFinite(P0) ? P0 * Pe_div_P0 : NaN;

  let Rspec = NaN;
  if (Number.isFinite(molarMass)) {
    Rspec = R_UNIVERSAL / molarMass;
  }
  const Ve =
    Number.isFinite(Rspec) && Number.isFinite(Te)
      ? Me * Math.sqrt(gamma * Rspec * Te)
      : NaN;

  let Cf = NaN;
  const Pe_div_P0_term = Math.pow(Pe_div_P0, (gamma - 1) / gamma);
  if (Number.isFinite(P0)) {
    const Pa_for_Cf = Number.isFinite(Pa) ? Pa : 0;
    const term1 = Math.sqrt(
      (2 * gamma * gamma) /
        (gamma - 1) *
        Math.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1)) *
        (1 - Pe_div_P0_term)
    );
    const term2 = ((Pe - Pa_for_Cf) / P0) * expansionRatio;
    Cf = term1 + term2;
  }

  let mdot = NaN;
  if (Number.isFinite(At) && Number.isFinite(Rspec) && Number.isFinite(T0) && Number.isFinite(P0)) {
    const chokeFactor = Math.pow(
      2 / (gamma + 1),
      (gamma + 1) / (2 * (gamma - 1))
    );
    mdot = P0 * At * Math.sqrt(gamma / (Rspec * T0)) * chokeFactor;
  }

  let thrust = NaN,
    thrustVac = NaN;
  if (Number.isFinite(Cf) && Number.isFinite(P0) && Number.isFinite(At)) {
    thrust = Cf * P0 * At;
    const term1_v = Math.sqrt(
      (2 * gamma * gamma) /
        (gamma - 1) *
        Math.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1)) *
        (1 - Pe_div_P0_term)
    );
    const term2_v = (Pe / P0) * expansionRatio;
    const Cf_vac = term1_v + term2_v;
    thrustVac = Cf_vac * P0 * At;
  }

  const g0 = 9.80665;
  let Isp = NaN;
  if (Number.isFinite(thrust) && Number.isFinite(mdot) && mdot > 0) {
    Isp = thrust / (mdot * g0);
  }

  document.getElementById("vExit").textContent = fmt(Ve);
  document.getElementById("pExit").textContent = fmt(Pe);
  document.getElementById("tExit").textContent = fmt(Te);
  document.getElementById("machExit").textContent = fmt(Me, 4);
  document.getElementById("pRatio").textContent = fmt(Pe_div_P0, 5);
  document.getElementById("tRatio").textContent = fmt(Te_div_T0, 5);
  document.getElementById("cf").textContent = fmt(Cf, 4);
  document.getElementById("thrust").textContent = fmt(thrust);
  document.getElementById("thrustVac").textContent = fmt(thrustVac);
  document.getElementById("isp").textContent = fmt(Isp, 3);
  document.getElementById("mdot").textContent = fmt(mdot);
}

document.getElementById("compute").addEventListener("click", calculateIsPerf);

// Get the hamburger button and the navigation menu elements
const menuToggle = document.getElementById('menuToggle');
const mainNav = document.getElementById('mainNav');

// Add a click event listener to the hamburger button
menuToggle.addEventListener('click', () => {
    // Toggle the 'mobile-open' class on the nav menu
    // This class will control the menu's visibility via CSS
    mainNav.classList.toggle('mobile-open');
    // Optionally, you can also toggle an 'active' class on the button itself for styling
    menuToggle.classList.toggle('active');
});   

// Close the menu when a link is clicked (for better UX on mobile)
mainNav.querySelectorAll('a').forEach(link => {
    link.addEventListener('click', () => {
        mainNav.classList.remove('mobile-open');
        menuToggle.classList.remove('active');
    });
}); 
// ---------- End of script.js ---------- 


