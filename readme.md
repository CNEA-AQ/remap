# remap

> Herramienta de remapeo e interpolación de arrays para fortran, basado en SCRIP. Pensado como módulo para modelos numéricos ó como herramienta para regrillar salidas de WRF.


## Dependencias:

- NetCDF Library
- [PROJ Library](https://proj.org/) 


## Remap methods:

Interpolation methods: (from coarser to finer grids)

| Method      | field-type  | Implemented? | 
|-------------|-------------|--------------| 
|Bilinear     | float       | [x]          | 
|Bicubic      | float       | [x]          | 
|Cubic-spline | float       | [ ]          | 
|Conserv. 1   | float       | [x]          | 
|Conserv. 2   | float       | [x]          | 
|Nearest N.   | integer     | [x]          | 

Remaping methods (from finer to coarser):

| Method      | field-type  | Implemented? | 
|-------------|-------------|--------------| 
| Average     | float       | [ ]          | 
| Weighted Avg| float       | [x]          | 
| Mode        | integer     | [ ]          | 



## Compilación

Ir al directorio `src/`:
```shell
cd src

```

Editar el archivo `Makefile` para indicar compilador, y la ubicación de la libreria netCDF y PROJ, y luego ejecutar el makefile:

```shell
make copy-links
make
``


## Ejecución

Cómo módulo, es necesario sólo llamarlo desde tu programa con:
```fortran
use SCRIP_mod 
```

y luego va a estár disponible la subrutina:

```fortran
call remap_field(arr1,arr2,g1,g2,method) 
```

donde:
  - arr1 es el array fuente (rango-1 ó rango-2) a ser remapeado.
  - arr2 es el array destino (rango-1 ó rango-2).
  - g1 es una estructura (*type*) que contiene la información de la grilla fuente.
  - g2 es una estructura (*type*) que contiene la información de la grilla destino.
  - méthod es el método de remapeo (bilinear, bicubic, diswgt, ó conservative).

## Próximas mejoras

- [ ] Terminar la versión `standalone`
- [ ] Dar soporte para remapeo de archivos de wrf directamente.
- [ ] Dar soporte a grillas logicamente no regulares.
- [ ] Producir los `grid_types` a partir de leer el archivo de entrada (wrf, cmaq o griddesc).

- [ ] Test rempaeo de WRF:
  + [x] a. agregar soporte a wildcards para iFile (<date>,<time>). 
  + [ ] b. calcular eta values en base a info de `ptop` y `e_levels`.
  + [ ] c.  
