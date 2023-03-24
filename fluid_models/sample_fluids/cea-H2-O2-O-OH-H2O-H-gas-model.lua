model = "CEAGas"

CEAGas = {
  mixtureName = 'h2-h-o2-o-oh-h2o',
  speciesList = {"H2", "O2", "O", "H", "OH", "H2O"},
  reactants = {H2=4/7, O2=3/7, O=0.0, H=0.0, OH=0.0, H2O=0.0},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}