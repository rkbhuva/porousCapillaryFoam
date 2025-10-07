# porousCapillaryFoam  

An extended OpenFOAM solver based on `twoPhaseEulerFoam` for simulating **capillary-driven multiphase transport in porous media** under isothermal conditions

- Developed and Tested with OpenFOAM v2306
## Features  
- **Capillary pressure models**  
  - Brooks–Corey  
  - van Genuchten  
  - Leverett function  

- **Relative permeability models**  
  - Brooks–Corey  
  - van Genuchten  

- **Physics extensions**  
  - Capillary pressure source term  
  - Darcy-type resistance for porous flow  
  - Isothermal formulation (energy equation disabled)

- This framework allows investigation of capillary-induced transport phenomena in multiphase porous media and provides a basis for further coupling with electrochemical or thermal models.

