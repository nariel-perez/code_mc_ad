# DIAGRAMA DE FLUJO - modulos_f90
## Simulación de Monte Carlo para Adsorción

```
═══════════════════════════════════════════════════════════════════════════════
                          ESTRUCTURA DE MÓDULOS
═══════════════════════════════════════════════════════════════════════════════

┌─────────────────────────────────────────────────────────────────────────────┐
│                           MÓDULOS PRINCIPALES                                │
└─────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────┐
│   InputParams        │  ◄─── Lee input.txt
│   (módulo)           │       Expone: P, dp, sigmetano, eps, ACEL, etc.
└──────────────────────┘
         │
         │ use
         ▼
┌──────────────────────┐
│   AdsorbateInput     │  ◄─── Lee MOLEC.DAT y archivos de moléculas
│   (módulo)           │       Expone: NMOLEC, maxAtoms, RX0, RY0, RZ0, etc.
└──────────────────────┘
         │
         │ use
         ▼
┌──────────────────────┐
│   SimulationData    │  ◄─── Variables globales de simulación
│   (módulo)           │       RX, RY, RZ, N, LOCATE, CNF, etc.
└──────────────────────┘

┌──────────────────────┐
│ PhysicalConstants    │  ◄─── Constantes físicas y eléctricas
│   (módulo)           │       ComputeConstants()
└──────────────────────┘

┌──────────────────────┐
│  RotationModule      │  ◄─── Tablas de rotación precomputadas
│   (módulo)           │       InitRotationTables(), GetRotationMatrix()
└──────────────────────┘

┌──────────────────────┐
│ EstructuraModule     │  ◄─── Estructura de la superficie
│   (módulo)           │       estructura()
└──────────────────────┘


═══════════════════════════════════════════════════════════════════════════════
                        FLUJO PRINCIPAL: program main
═══════════════════════════════════════════════════════════════════════════════

START
  │
  ├─► call cpu_time(starttime)
  │
  ├─► [LECTURA DE ENTRADA]
  │   │
  │   ├─► call read_input('input.txt')          [InputParams]
  │   │   └─► Lee: P, dp, sigmetano, eps, ACEL, diel, BCX/BCY/BCZ,
  │   │           T, mat, nam, isot, ijpasos, ikpasos, mult2,
  │   │           ensemble, NCELLMAT, NESTADO, canonicalmolecules
  │   │
  │   └─► call print_params()                   [InputParams]
  │       └─► Imprime parámetros y valida NCELLMAT <= 50
  │
  ├─► [PREPARACIÓN INICIAL]
  │   │
  │   ├─► if (ensemble == 3) then
  │   │   └─► call namd1()                       [¿namd1.f90?]
  │   │       └─► ensemble = 1
  │   │
  │   ├─► allocate(uads(...))                   [SimulationData]
  │   │
  │   ├─► [Construcción de presiones logarítmicas]
  │   │   └─► p_vals(1..isot+1) = valores logarítmicos
  │   │
  │   ├─► [Unidades reducidas]
  │   │   └─► SIGMA = sigmetano/acel
  │   │       TEMP = T/eps
  │   │       PRED = P*SIGMA³/eps
  │   │       RCUT = 10*SIGMA
  │   │       VOL = XMAX*YMAX*ZMAX
  │   │
  │   ├─► call computeconstants(...)             [PhysicalConstants]
  │   │   └─► Calcula: AK, AE0, FACTORELEC, FCLEC, etc.
  │   │
  │   └─► call initrotationtables()              [RotationModule]
  │       └─► Precomputa TABLE_SEN y TABLE_COS
  │
  ├─► [LECTURA DE ADSORBATOS]
  │   │
  │   └─► call read_adsorbates('MOLEC.DAT', sigma, 1.0e-7)  [AdsorbateInput]
  │       │
  │       ├─► Lee MOLEC.DAT → NMOLEC
  │       ├─► Para cada molécula:
  │       │   ├─► Lee nombre del archivo
  │       │   ├─► Lee: nAtomFile, fracción X, NMIN2, NMAXI2
  │       │   └─► Determina maxAtoms
  │       │
  │       ├─► allocate(X, NATOM, NMIN2, NMAXI2, NTOTAL, ...)
  │       ├─► allocate(RX0, RY0, RZ0, NATOMKIND, NSYM)
  │       │
  │       └─► Para cada molécula:
  │           ├─► Lee coordenadas (x1, y1, z1)
  │           ├─► Escala: RX0 = x1/sigmetano * sigma
  │           └─► Guarda NATOMKIND y NSYM
  │
  ├─► [ALOCACIÓN DE ARRAYS PRINCIPALES]
  │   │
  │   ├─► allocate(Z(NMOLEC))
  │   ├─► allocate(N(NMOLEC))                    [SimulationData]
  │   ├─► allocate(EPSI, SIGM, Q(maxAtoms))
  │   ├─► allocate(RX1, RY1, RZ1(maxAtoms))
  │   ├─► allocate(RX, RY, RZ(5000, maxAtoms, NMOLEC))
  │   ├─► allocate(USS(5000, maxAtoms, maxAtoms))
  │   └─► allocate(CNF(...))                      [SimulationData]
  │
  ├─► [INICIALIZACIÓN DE SUPERFICIE Y POTENCIALES]
  │   │
  │   ├─► open('SALIDAACTIVADO-100.TXT')
  │   ├─► open('PERFILES.TXT')
  │   │
  │   ├─► call estructura(eps, nam, sigma, sigmetano, NC, diel)  [EstructuraModule]
  │   │   └─► Lee archivo de superficie 'nam'
  │   │       Calcula estructura cristalina
  │   │
  │   ├─► call potencialff(...)                  [PotencialFF.f90]
  │   │   └─► Construye tablas de potencial fuerza-fuerza
  │   │
  │   └─► call potencial(...)                    [Potencial.f90]
  │       └─► Construye tablas de potencial adsorbato-superficie
  │
  ├─► [INICIALIZACIÓN SEGÚN ENSEMBLE]
  │   │
  │   ├─► if (ensemble == 2) then                [Gran Canónico - arranque en cero]
  │   │   └─► N(:) = 0, NTOTAL(:) = 0
  │   │       LOCATE(:,:) = 0
  │   │       V = 0, VG = 0, VA = 0
  │   │
  │   └─► if (ensemble < 2) then                 [Canónico o GC con restart]
  │       │
  │       ├─► open('initconf.txt')
  │       ├─► read V, VG, VA
  │       ├─► read nmolec2
  │       │
  │       └─► Para cada tipo de molécula:
  │           ├─► read natom2, ncantmol
  │           └─► Para cada molécula:
  │               ├─► read coordenadas → RX1, RY1, RZ1
  │               └─► call add(molkind)          [Add.f90]
  │                   └─► Actualiza LOCATE y arrays RX/RY/RZ
  │       │
  │       └─► [Cálculo de energía inicial]
  │           └─► Para cada molécula:
  │               ├─► ipull = locate(J,I)
  │               ├─► call adpotout(ipull, I, deltva)  [Adpotout.f90]
  │               └─► call potout(ipull, I, deltv)     [Potout.f90]
  │
  ├─► [BUCLE PRINCIPAL: ISOTERMAS]
  │   │
  │   └─► do IPASOS = 1, isot                    [isot puntos de isoterma]
  │       │
  │       ├─► [Limpieza de estadísticas]
  │       │   └─► CNF(:,:,:,:,:) = 0
  │       │
  │       ├─► open('CONFIG###.xyz', 'CONFIG###.TXT')
  │       │
  │       ├─► [Cálculo de actividades Z]
  │       │   └─► if (NESTADO == 1) then
  │       │           Z(INMOLEC) = X*P/(8.3144*T)*6.023E-7*auvol
  │       │       else
  │       │           Z(INMOLEC) = X*P*6.023E-4*(ACEL³)*VOL
  │       │
  │       ├─► [Inicialización de acumuladores]
  │       │   └─► U = 0, UG = 0, UA = 0, UN = 0, etc.
  │       │
  │       ├─► [BUCLE DE PROMEDIO: JPASOS]
  │       │   │
  │       │   └─► do JPASOS = 1, ijpasos         [ijpasos pasos para promediar]
  │       │       │
  │       │       ├─► V = V/eps, VG = VG/eps, VA = VA/eps  [Unidades reducidas]
  │       │       ├─► MULT = mult2 (si JPASOS == 1)
  │       │       │
  │       │       ├─► [BUCLE DE MOVIMIENTOS: KPASOS]
  │       │       │   │
  │       │       │   └─► do KPASOS = 1, ikpasos*MULT
  │       │       │       │
  │       │       │       ├─► IJ = ranf(DUMMY)*3 + 1        [Ranf.f90]
  │       │       │       │
  │       │       │       ├─► if (ensemble == 0) goto 30    [Solo MOVE]
  │       │       │       │
  │       │       │       └─► goto (10, 20, 30, 35) IJ
  │       │       │           │
  │       │       │           ├─► [10] call in(...)         [In.f90]
  │       │       │           │   │
  │       │       │           │   ├─► Elige tipo de molécula (MOLKIND)
  │       │       │           │   ├─► Genera posición aleatoria (RXNEW, RYNEW, RZNEW)
  │       │       │           │   ├─► Carga coordenadas base: RX0 → RXBE
  │       │       │           │   ├─► call GetRotationMatrix(...)  [RotationModule]
  │       │       │           │   │   └─► Calcula matriz de rotación R(3,3)
  │       │       │           │   ├─► Rota molécula: RX1 = R * RXBE + RXNEW
  │       │       │           │   ├─► call potin(...)        [Potin.f90]
  │       │       │           │   │   └─► Calcula DELTV (energía adsorbato-adsorbato)
  │       │       │           │   ├─► call adpotin(...)     [Adpotin.f90]
  │       │       │           │   │   └─► Calcula DELTVA (energía adsorbato-superficie)
  │       │       │           │   ├─► Criterio de aceptación Metropolis
  │       │       │           │   └─► if (aceptado) then
  │       │       │           │       └─► call add(molkind)  [Add.f90]
  │       │       │           │           └─► Actualiza LOCATE, RX/RY/RZ, N
  │       │       │           │
  │       │       │           ├─► [20] call out(...)        [Out.f90]
  │       │       │           │   │
  │       │       │           │   ├─► Elige tipo de molécula (MOLKIND)
  │       │       │           │   ├─► Elige molécula aleatoria: NLOC → IPULL
  │       │       │           │   ├─► call potout(IPULL, MOLKIND, DELTV)  [Potout.f90]
  │       │       │           │   │   └─► Calcula DELTV al eliminar
  │       │       │           │   ├─► call adpotout(IPULL, MOLKIND, DELTVA)  [Adpotout.f90]
  │       │       │           │   │   └─► Calcula DELTVA al eliminar
  │       │       │           │   ├─► Criterio de aceptación Metropolis
  │       │       │           │   └─► if (aceptado) then
  │       │       │           │       └─► call remove(NLOC, IPULL, MOLKIND)  [Remove.f90]
  │       │       │           │           └─► Actualiza LOCATE, N
  │       │       │           │
  │       │       │           ├─► [30] call move(...)        [Move.f90]
  │       │       │           │   │
  │       │       │           │   ├─► Elige molécula aleatoria: NLOC → IPULL
  │       │       │           │   ├─► Guarda posición antigua
  │       │       │           │   ├─► Genera desplazamiento aleatorio
  │       │       │           │   ├─► call GetRotationMatrix(...)  [RotationModule]
  │       │       │           │   ├─► Rota y desplaza molécula
  │       │       │           │   ├─► call potout(...) + potin(...)  [Potout/Potin]
  │       │       │           │   ├─► call adpotout(...) + adpotin(...)  [Adpotout/Adpotin]
  │       │       │           │   ├─► Criterio de aceptación Metropolis
  │       │       │           │   └─► if (aceptado) then
  │       │       │           │       └─► Actualiza RX/RY/RZ
  │       │       │           │
  │       │       │           └─► [35] call change(...)     [change.f90]
  │       │       │               └─► Cambio de orientación/conformación
  │       │       │
  │       │       ├─► [Después del bucle KPASOS]
  │       │       │   │
  │       │       │   ├─► V = V*eps, VG = VG*eps, VA = VA*eps  [Volver a unidades físicas]
  │       │       │   │
  │       │       │   ├─► Acumulación de estadísticas
  │       │       │   │   └─► NTOTAL(I) = NTOTAL(I) + N(I)
  │       │       │   │
  │       │       │   └─► call estadistica(NCELLMAT)  [estadistica.f90]
  │       │       │       │
  │       │       │       └─► Para cada molécula y átomo:
  │       │       │           ├─► ICNF = INT(RX * NCELLMAT)
  │       │       │           ├─► JCNF = INT(RY * NCELLMAT)
  │       │       │           ├─► KCNF = INT(RZ * NCELLMAT)
  │       │       │           └─► CNF(ICNF, JCNF, KCNF, I, NATOMKINDI) += 1
  │       │       │
  │       │       └─► Acumulación de energías
  │       │           └─► U = U + V, UG = UG + VG, UA = UA + VA
  │       │               UN = UN + V*NTOTALB, etc.
  │       │
  │       ├─► [Después del bucle JPASOS]
  │       │   │
  │       │   ├─► Guardar configuración en 'initconf.txt'
  │       │   │   └─► write V, VG, VA
  │       │   │       write coordenadas de todas las moléculas
  │       │   │
  │       │   ├─► Cálculo de promedios
  │       │   │   └─► ANPROM(I) = NTOTAL(I) / ijpasos
  │       │   │
  │       │   ├─► Cálculo de calores isostéricos
  │       │   │   └─► CALOR = 8.3144*T - ((UN1 - U1*AN1)/ANN)*8.31
  │       │   │       CALORG = (UNG1 - UG1*AN1)/ANN*8.31
  │       │   │       CALORA = -((UNA1 - UA1*AN1)/ANN)*8.31
  │       │   │
  │       │   ├─► Escritura de configuración XYZ
  │       │   │   └─► write('CONFIG###.xyz')
  │       │   │       └─► Escribe superficie + moléculas adsorbidas
  │       │   │
  │       │   └─► Escritura de estadísticas espaciales CNF
  │       │       └─► write('CNF###-I-NATOMKINDI.TXT')
  │       │
  │       └─► Actualizar presión: P = p_vals(IPASOS + 1)
  │
  ├─► [FINALIZACIÓN]
  │   │
  │   ├─► close archivos de salida
  │   ├─► deallocate arrays
  │   └─► call cpu_time(endtime)
  │       └─► write tiempo de ejecución
  │
  END


═══════════════════════════════════════════════════════════════════════════════
                        DEPENDENCIAS ENTRE MÓDULOS
═══════════════════════════════════════════════════════════════════════════════

Main.f90
  ├─► use InputParams          → read_input(), print_params()
  ├─► use SimulationData       → Variables globales (RX, RY, RZ, N, LOCATE, CNF, ...)
  ├─► use AdsorbateInput       → read_adsorbates(), NMOLEC, maxAtoms, RX0, RY0, RZ0, ...
  ├─► use PhysicalConstants    → ComputeConstants()
  ├─► use EstructuraModule     → estructura()
  └─► use RotationModule       → InitRotationTables(), GetRotationMatrix()

InputParams.f90
  └─► (módulo independiente)

AdsorbateInput.f90
  └─► use InputParams          → solo: sigmetano

SimulationData.f90
  └─► (módulo independiente)

PhysicalConstants.f90
  └─► (módulo independiente)

RotationModule.f90
  └─► (módulo independiente)

EstructuraModule.f90
  └─► (probablemente usa SimulationData y PhysicalConstants)


═══════════════════════════════════════════════════════════════════════════════
                        SUBRUTINAS DE MOVIMIENTO MC
═══════════════════════════════════════════════════════════════════════════════

In.f90 (CREACIÓN)
  ├─► use InputParams          → acel, acelx, acely, acelz, bcx, bcy, bcz
  ├─► use AdsorbateInput       → RX0, RY0, RZ0, NMOLEC, NATOM
  ├─► use SimulationData      → RX1, RY1, RZ1, EXNEW, EYNEW, EZNEW, ANX, ANGY, ANZ, N
  ├─► use RotationModule       → GetRotationMatrix()
  ├─► call GetRotationMatrix() → Calcula matriz de rotación
  ├─► call potin()            → [Potin.f90] Energía adsorbato-adsorbato
  ├─► call adpotin()          → [Adpotin.f90] Energía adsorbato-superficie
  └─► call add()              → [Add.f90] Inserta molécula

Out.f90 (DESTRUCCIÓN)
  ├─► use SimulationData      → N, LOCATE
  ├─► use AdsorbateInput       → NMOLEC
  ├─► call potout()           → [Potout.f90] Energía adsorbato-adsorbato
  ├─► call adpotout()         → [Adpotout.f90] Energía adsorbato-superficie
  └─► call remove()           → [Remove.f90] Elimina molécula

Move.f90 (DESPLAZAMIENTO)
  ├─► use InputParams          → acel, acelx, acely, acelz, bcx, bcy, bcz
  ├─► use AdsorbateInput       → RX0, RY0, RZ0, NATOM, NMOLEC
  ├─► use SimulationData      → RX, RY, RZ, RX1, RY1, RZ1, N, LOCATE, ANX, ANGY, ANZ, EXNEW, EYNEW, EZNEW
  ├─► use RotationModule       → GetRotationMatrix()
  ├─► call GetRotationMatrix() → Calcula matriz de rotación
  ├─► call potout() + potin()  → [Potout/Potin] Energía adsorbato-adsorbato
  └─► call adpotout() + adpotin() → [Adpotout/Adpotin] Energía adsorbato-superficie

change.f90 (CAMBIO DE ORIENTACIÓN/CONFORMACIÓN)
  └─► (similar a Move pero para cambios conformacionales)


═══════════════════════════════════════════════════════════════════════════════
                        ARCHIVOS DE ENTRADA/SALIDA
═══════════════════════════════════════════════════════════════════════════════

ENTRADA:
  ├─► input.txt               → [InputParams.read_input()]
  │   └─► P, dp, sigmetano, eps, ACEL, acelx, acely, acelz, diel,
  │       BCX, BCY, BCZ, T, mat, nam, isot, ijpasos, ikpasos,
  │       mult2, ensemble, NCELLMAT, NESTADO, canonicalmolecules
  │
  ├─► MOLEC.DAT               → [AdsorbateInput.read_adsorbates()]
  │   └─► NMOLEC
  │       └─► Para cada molécula: nombre del archivo
  │
  ├─► [archivos de moléculas]  → [AdsorbateInput.read_adsorbates()]
  │   └─► nAtomFile, fracción X, NMIN2, NMAXI2
  │       └─► Para cada átomo: x, y, z, ikind, ns
  │
  ├─► [archivo de superficie]  → [estructura()]
  │   └─► nombre dado por variable 'nam'
  │
  └─► initconf.txt            → [Main - solo si ensemble < 2]
      └─► V, VG, VA
          └─► nmolec2
              └─► Para cada tipo: natom2, ncantmol, coordenadas

SALIDA:
  ├─► SALIDAACTIVADO-100.TXT   → [Main] Resultados principales de isoterma
  ├─► PERFILES.TXT            → [Main] Perfiles de densidad
  ├─► CONFIG###.xyz           → [Main] Configuraciones en formato XYZ
  ├─► CONFIG###.TXT           → [Main] Configuraciones en formato texto
  ├─► CNF###-I-NATOMKINDI.TXT → [Main] Estadísticas espaciales (CNF)
  ├─► Ener###.TXT             → [Main] Energías durante la simulación
  └─► initconf.txt            → [Main] Configuración final (para restart)


═══════════════════════════════════════════════════════════════════════════════
                        LEYENDA DE SÍMBOLOS
═══════════════════════════════════════════════════════════════════════════════

  ──►  Flujo de ejecución
  ├─►  Rama del flujo
  └─►  Final de rama
  │    Continuación vertical
  ◄─── Indica origen de datos/archivo
  [ ]  Nombre de módulo o archivo
  ()   Parámetros o argumentos
  call Función/subrutina llamada
