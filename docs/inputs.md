
# Input Options

## The Universe Class

The `Universe` class provides information about the background universe to the rest of the `DeepPhase` code.

### Parameters

| Parameter | Type   | Description | Default |
|-----------|--------|-------------|---------|
| `T0`      | double | Current temperature of the universe in GeV | 2.41e-13 |
| `H0`      | double | Hubble rate today (in GeV) | 1.45e-42 |
| `g0`      | double | Relativistic degrees of freedom today | 3.91 |
| `Ts`      | double | Temperature during phase transition (GeV) | 100.0 |
| `Hs`      | double | Hubble rate at phase transition (GeV) | 1.41e-14 |
| `gs`      | double | Relativistic degrees of freedom at transition | 106.75 |

### Constructors

```cpp
Universe myUniverse; // Default constructor
Universe myUniverse(100.0, 1.41e-14, 106.75); // Ts, Hs, gs (current-time values set to defaults)
Universe myUniverse(2.41e-13, 100.0, 1.45e-42, 1.41e-14, 3.91, 106.75); // T0, Ts, H0, Hs, g0, gs
```

---

## The PTParams Class

The `PTParams` class provides information from the cosmological phase transition to the rest of the `DeepPhase` code.

### Parameters

| Parameter    | Type        | Description | Default |
|--------------|-------------|-------------|---------|
| `vw`         | double      | Bubble wall velocity | 0.8 |
| `alN`        | double      | Strength parameter $\alpha$ at nucleation | 0.1 |
| `betaH`      | double      | Normalized transition rate $\beta/H$ | 10.0 |
| `dtauH`      | double      | Normalized sound wave lifetime | 1.0 |
| `wNeN_rat`   | double      | Enthalpy/energy density ratio | 4./3. |
| `nuc_type`   | std::string | Nucleation type ("exp") | "exp" |
| `un`         | Universe    | Universe object | Universe() |

### Constructors

```cpp
PTParams myParams; // Default constructor
PTParams myParams(0.8, 0.1, 10.0, 1.0, 4./3., "exp", myUniverse); // All parameters
```
