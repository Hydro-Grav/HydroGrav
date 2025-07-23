# Fluid Profile Inputs
* `vw` : The bubble wall velocity, with the default value `0.8`.
* `alN` : The $\alpha$ strength parameter of the phase transition, with default value `alN = 0.1`.
* `beta` : The inverse time scale $\beta$ of the phase transition in GeV, with default value `1.0`.
* `dtau` : The duration of the phase transition in GeV$^{-1}$, with default value $10.0$. Typically taken to be $\beta^{-1}$.
* `wNeN_rat` : The enthalpy density to energy density ratio of the phase transition, with default value `4./3.` determined from the bag model.
* `nuc_type` : The nucleation type, which currently only supports the option `"exp"`.

# Universe Inputs
* `T0` : Current temperature of the universe in GeV. The default value of `T0 = 2.41e-13` corresponds to $T_0 = 2.725$ K.
* `H0` : Hubble rate measured today. The default value of `H0 = 1.45e-42` corresponds to the known value $H_0 = 100h$ km/s/Mpc, where we take $h = 0.678$ (Planck 2018).
* `g0` : Current relativistic degrees of freedom, with the default value `g0 = 3.91`.
* `Ts` : Temperature of the universe during the phase transition in GeV, taken to be either the nucleation or percolation temperatures of the transition. The default value `Ts = 100.0` corresponds to the EW scale.
* `Hs` : Hubble rate at the time of the phase transition. The default value of `Hs = 1.41e-14` is given by $H_s^2 = 8 \pi G \rho_s /3 = 8 \pi^3 G g_s T_s^4/90$.
* `gs` : Relativistic degrees of freedom during the phase transition. The default value of `gs = 106.75` corresponds to the SM at `Ts`.