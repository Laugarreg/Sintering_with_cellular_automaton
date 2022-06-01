MODULE variables
implicit none
! *********************************************
! *********************************************
! 			DECLARACION DE VARIABLES
! *********************************************
! *********************************************
INTEGER :: n_mat		! Tamanio de la caja 
REAL(8) :: rho			! Densidad teorica
REAL(8) :: radio		! Radio de los circulos
INTEGER :: n_int		! Numero de intentos para colocar circulos

INTEGER, DIMENSION(:,:), ALLOCATABLE :: matriz	! Matriz del sistema
REAL(8), DIMENSION(:,:), ALLOCATABLE :: circulos
	! Coordenas del centro de los circulos y radio
	
INTEGER :: circ_colocados	! Numero de circulos colocados
INTEGER :: px_materia		! Numero de celdas con materia en matriz
REAL(8) :: densidad_real	! Densidad real en cada iteracion
INTEGER :: intentos			! Numero de intentos en cada iteracion
REAL(8) :: rnd				! Auxiliar para numero aleatorio [0,1)
REAL(8) :: x_new, y_new		! Posicion para los centros de los circulos
REAL(8) :: x_rel, y_rel		! Posicion relativa
REAL(8) :: dist				! Distancia entre centros de dos circulos

!  FIN DEL MODULO
END MODULE variables


PROGRAM CIRCULOS_CON_SOLAPAMIENTO
USE variables

! DECLARACION DE VARIABLES LOCALES
implicit none
INTEGER :: i		! Indice para bucles DO


! *********************************************
! *********************************************
! 			INPUTS / INICIALIZACION
! *********************************************
! *********************************************
! Lectura de inputs desde fichero externo
WRITE(*,*) 'Leyendo inputs'
OPEN(UNIT=11, FILE='circulos-con-solapamiento.in', STATUS='OLD')
	READ(11,*) n_mat
	READ(11,*) rho
	READ(11,*) radio
	READ(11,*) n_int
CLOSE(11)

WRITE(*,*) 'Tamanio del sistema: ', n_mat
WRITE(*,*) 'Densidad teorica: ', rho
WRITE(*,*) 'Radio de los circulos: ', radio
WRITE(*,*) 'Numero de intentos para colocar circulos: ', n_int
WRITE(*,*) ' '

! Inicializacion de magnitudes
ALLOCATE(matriz(0:n_mat-1, 0:n_mat-1))
matriz  = 0
WRITE(*,*) 'Matriz inicializada'

ALLOCATE(circulos(1:n_mat*n_mat, 3))
circulos = 0.0d0
WRITE(*,*) 'Circulos inicializados'

circ_colocados = 0
px_materia = 0
densidad_real = 0.0d0
intentos = 0
WRITE(*,*) 'Variables auxiliares inicializadas'
WRITE(*,*) ' '


! *********************************************
! *********************************************
! 				CODIGO PRINCIPAL
! *********************************************
! *********************************************
! INICIALIZACION DE FICHEROS DE ESCRITURA
WRITE(*,*) 'Inicializacion de ficheros de escritura'
OPEN(UNIT=21, FILE='densidad-final.dat', STATUS='UNKNOWN')
OPEN(UNIT=22, FILE='coordenadas-circulos.dat', STATUS='UNKNOWN')
OPEN(UNIT=23, FILE='matriz-sistema.dat', STATUS='UNKNOWN')


! CALCULO PRINCIPAL
WRITE(*,*) 'Comienzo del calculo principal'

! Se coloca el primer circulo
CALL RANDOM_NUMBER(rnd)
x_new = DBLE(n_mat) * rnd
CALL RANDOM_NUMBER(rnd)
y_new = DBLE(n_mat) * rnd

circ_colocados = circ_colocados + 1
circulos(circ_colocados,1) = x_new
circulos(circ_colocados,2) = y_new
circulos(circ_colocados,3) = radio
CALL UPDATE_MATRIX

! Se intentan colocar circulos hasta alcanzar la densidad teorica o 
! el maximo numero de intentos permitidos
DO WHILE (densidad_real.LT.rho .AND. intentos.LT.n_int)
	
	! Se genera una posicion para el centro del nuevo circulo
	CALL RANDOM_NUMBER(rnd)
	x_new = DBLE(n_mat) * rnd
	CALL RANDOM_NUMBER(rnd)
	y_new = DBLE(n_mat) * rnd
	
	! Se añade el circulo, se actualiza la matriz y se  recalcula la 
	! densidad_real
	intentos = 0
	circ_colocados = circ_colocados + 1
	circulos(circ_colocados,1) = x_new
	circulos(circ_colocados,2) = y_new
	circulos(circ_colocados,3) = radio
	
	CALL UPDATE_MATRIX
	WRITE(*,*) 'Densidad: ', densidad_real
	
ENDDO
WRITE(*,*) 'Fin del calculo principal'
WRITE(*,*) ' '


! ESCRITURA DE RESULTADOS EN FICHERO
! Densidad final
WRITE(21,*) densidad_real

! Centros y radios de los circulos
WRITE(22,*) circ_colocados
WRITE(22,*) 'Circulos colocados sin solapar'
DO i=1,circ_colocados
	WRITE(22,*) circulos(i,:)
ENDDO

! Matriz del sistema
DO i=0,n_mat-1
	WRITE(23,*) matriz(i,:)
ENDDO
WRITE(*,*) 'Resultados salvados'


! CERRADO DE FICHEROS DE ESCRITURA
WRITE(*,*) 'Cerrado de ficheros de escritura'
CLOSE(21)
CLOSE(22)
CLOSE(23)
WRITE(*,*) ' '


! VISUALIZACION
CALL system('gnuplot -p volcado.plt')

WRITE(*,*) '---------------'
WRITE(*,*) 'Done!'


CONTAINS
! *********************************************
! 		SUBROUTINE: UPDATE_MATRIX
! *********************************************
SUBROUTINE UPDATE_MATRIX
USE variables, ONLY : n_mat, x_new, y_new, radio, matriz, &
						px_materia, densidad_real

! DECLARACION DE VARIABLES LOCALES
implicit none
INTEGER :: i, j 			! Indices para bucles DO
REAL(8) :: x_celda, y_celda	! Centro de cada celda de la matriz
REAL(8) :: x_rel_aux, y_rel_aux		! Posiciones relativas auxiliares
REAL(8) :: dist_aux	! Distancia auxiliar para la ocupacion de celdas

! SUBRUTINA
! Se comprueba que pixeles de la matriz se convierten en materia
DO i=0,n_mat-1
	DO j=0,n_mat-1
	
		! Calculo del centro de cada celda
		x_celda = DBLE(i) + 0.5
		y_celda = DBLE(j) + 0.5
		
		! Posiciones relativas
		x_rel_aux = x_celda - x_new
		y_rel_aux = y_celda - y_new
		
		! Correcion PBC
		x_rel_aux = x_rel_aux - DBLE(n_mat)*DNINT(x_rel_aux/DBLE(n_mat)) 
		y_rel_aux = y_rel_aux - DBLE(n_mat)*DNINT(y_rel_aux/DBLE(n_mat)) 
		
		! Distancia del centro de la celda al centro del nuevo circulo
		dist_aux = SQRT(x_rel_aux**2.0d0 + y_rel_aux**2.0d0)
		IF (dist_aux.LT.radio .AND. matriz(i,j).NE.1) THEN
			matriz(i,j) = 1
			px_materia = px_materia + 1
		ENDIF
	ENDDO
ENDDO

! Se actualiza la densidad_real
densidad_real = DBLE(px_materia)/(n_mat*n_mat)

! FIN DE LA SUBRUTINA
END SUBROUTINE UPDATE_MATRIX

END PROGRAM CIRCULOS_CON_SOLAPAMIENTO
