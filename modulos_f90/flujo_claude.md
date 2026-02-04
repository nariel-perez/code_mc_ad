# Resumen del Flujo de Trabajo - Simulación GCMC

Simulación de **Monte Carlo Gran Canónico (GCMC)** para adsorción de gases en superficies porosas.

---

## 1. Inicialización (`Main.f90`)

```
input.txt → InputParams    (parámetros de simulación)
MOLEC.DAT → AdsorbateInput (estructura de moléculas adsorbato)
```

- Carga parámetros: presión, temperatura, dimensiones de celda, tipo de ensemble
- Calcula constantes físicas y eléctricas (`PhysicalConstants`)
- Precalcula tablas trigonométricas (`RotationModule`)

---

## 2. Construcción de la Superficie (`Estructura.f90`)

- Lee archivo de superficie (nombre en `nam`, p. ej. `carbo.txt`)
- Filtra átomos dentro de la caja de simulación
- Convierte coordenadas a unidades reducidas
- Escribe archivo truncado: **`<base>_truncado.txt`** (base = nombre sin extensión; ej. `carbo.txt` → `carbo_truncado.txt`). Main usa este archivo para escribir las configuraciones XYZ.

---

## 3. Tablas de Potencial (precálculo)

| Archivo | Propósito |
|---------|-----------|
| `PotencialFF.f90` | Potencial fluido-fluido → tabla `USS(dist, tipo1, tipo2)` |
| `Potencial.f90` | Potencial adsorbato-superficie → tabla `UADS(x,y,z,tipo)` |

---

## 4. Bucle Principal de Simulación

```
IPASOS (isotermas) → JPASOS (promedios) → KPASOS (movimientos MC)
```

En cada paso MC se elige aleatoriamente:

| Código | Subrutina | Acción |
|--------|-----------|--------|
| 10 | `In.f90` | **Crear** molécula (inserción) |
| 20 | `Out.f90` | **Destruir** molécula (eliminación) |
| 30 | `Move.f90` | **Mover** molécula (traslación/rotación) |
| 35 | `change.f90` | **Cambiar** tipo de molécula |

---

## 5. Criterio de Aceptación Metropolis

Cada operación calcula:
- `DELTV` → cambio de energía fluido-fluido (`Potin/Potout`)
- `DELTVA` → cambio de energía adsorbato-superficie (`Adpotin/Adpotout`)
- Acepta/rechaza según: `exp(-β×ΔE) > random`

---

## 6. Gestión de Moléculas

- `Add.f90`: Inserta molécula en arrays `RX`, `RY`, `RZ` y actualiza `LOCATE`
- `Remove.f90`: Compacta el array `LOCATE` tras eliminación

---

## 7. Estadísticas y Salida

- `estadistica.f90`: Acumula histogramas espaciales (`CNF`)
- Calcula calores isostéricos
- Escribe configuraciones (`.xyz`, `.TXT`)
- Guarda `initconf.txt` para restart

---

## Diagrama Simplificado

```
┌──────────────┐     ┌─────────────────┐     ┌──────────────┐
│  Entrada     │ ──► │  Precálculo     │ ──► │  Monte Carlo │
│  input.txt   │     │  Potenciales    │     │  GCMC Loop   │
│  MOLEC.DAT   │     │  USS, UADS      │     │  IN/OUT/MOVE │
│  superficie  │     │                 │     │              │
└──────────────┘     └─────────────────┘     └──────┬───────┘
                                                    │
                     ┌─────────────────┐            │
                     │  Salida         │ ◄──────────┘
                     │  isotermas      │
                     │  CONFIG*.xyz    │
                     │  CNF*.TXT       │
                     └─────────────────┘
```

---

## Notas

- El código simula adsorción de múltiples especies moleculares sobre materiales porosos (probablemente carbón activado), calculando isotermas de adsorción y calores isostéricos.
- **Compilación**: script `run2` (bash run2). Genera `simulation.exe` y elimina al final `.o` y `.mod`.
- Documentación detallada: `DIAGRAMA_FLUJO.md`, `MEJORAS.md`, `ESTADO_DEL_CODIGO.md`.
