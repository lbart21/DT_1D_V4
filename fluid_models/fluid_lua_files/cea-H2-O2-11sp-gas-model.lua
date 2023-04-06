model = "CEAGas"

CEAGas = {
  mixtureName = 'h2-o2-h-o-h2o-ho2-oh-h2o2-oh-ar-he-n2',
  speciesList = {"H2", "O2", "H", "O", "H2O", "HO2", "OH", "H2O2", "Ar", "He", "N2"},
  reactants = {H2=0.005, O2=0.005, N2=0.99},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}
