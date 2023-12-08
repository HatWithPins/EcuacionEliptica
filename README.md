# EcuacionEliptica

:es: \
Implementaci�n del m�todo de las diferencias finitas para resolver la ecuaci�n de Laplace de medio c�rculo en coordenadas polares y cartesianas. Usa CMake para la compilaci�n.

Para compilar en la l�nea de comandos:

`cmake` \
`cmake --build <carpeta-de-salida> --target EcuacionEliptica`

Si se usa un IDE tipo Visual Studio, abre la carpeta como un proyecto y compila como cualquier otro proyecto. Aseg�rate de que tienes la extensi�n para CMake.

Una vez se ha compilado, hay que pasarle las condiciones de contorno y el n�mero de particiones.

`EcuacionEliptica boundary_bottom=0 boundary_arc=1 N=3 M=6`

Si los argumentos tienen alg�n problema, devuelve error.

## Estructura de los archivos.
- EcuacionElpitica contiene el main y las comprobaciones de los argumentos.
- Solver contiene las funciones para plantear el problema y resolverlo.

:uk: \
Implementation of finite differences method for solving Laplace's equation for a half circle in polar and cartesian coordinates. Builds using cmake.

To build using command line:

`cmake` \
`cmake --build <output-folder> --target EcuacionEliptica`

In case you use and IDE like Visual Studio, open the folder as a project and build it like any other project. Ensure that your IDE has the extensions for CMake.

Once is built, it needs the boundary conditions and number of partitions.

`EcuacionEliptica boundary_bottom=0 boundary_arc=1 N=3 M=6`

If some of the arguments has any issue, program will return an error.

## Folder structure.
- EcuacionElpitica contains the main function and is in charge of checking parsed arguments.
- Solver contains the functions for possing and solving the problem.