model = "CEAGas"

CEAGas = {
  mixtureName = 'h2-o2-h2o',
  speciesList = {"H2", "O2", "H2O"},
  reactants = {H2=77.7/144.6, O2=66.9/144.6, H2O=0.0},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}