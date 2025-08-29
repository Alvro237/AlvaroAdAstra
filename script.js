function calculateLift() {
  // Get values from inputs
  const V = parseFloat(document.getElementById("velocity").value);
  const rho = parseFloat(document.getElementById("density").value);
  const S = parseFloat(document.getElementById("area").value);
  const Cl = parseFloat(document.getElementById("cl").value);

  // Check if inputs are valid
  if (isNaN(V) || isNaN(rho) || isNaN(S) || isNaN(Cl)) {
    alert("Please enter all values!");
    return;
  }

  // Lift formula: L = 0.5 * rho * V^2 * S * Cl
  const L = 0.5 * rho * V * V * S * Cl;

  // Show result
  document.getElementById("result").textContent = L.toFixed(2);
}
