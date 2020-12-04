      !-------------------------------------------------------------
      ! function kinetic_energy: Calculate kinetic energy based on velocities
      !
      ! Input:
      !   natoms: number of atoms in the system
      !   vel1: (natoms,3) velocities of all atoms
      !   massO: Mass of 1 Ar atom
      !
      ! Output:
      !   ke: Total kinetic energy in the system
      !-------------------------------------------------------------
      function kinetic_energy(nMol,nAtomPer,velold,vel1,mass) result(ke)
        IMPLICIT NONE

        INTEGER, intent(in):: nMol        ! number of molecule
        INTEGER, intent(in):: nAtomPer        ! number of atoms per mol
        REAL*8, intent(in) :: velold(nMol,nAtomPer,3)  ! velocity array
        REAL*8, intent(in) :: vel1(nMol,nAtomPer,3)  ! velocity array

        REAL*8, intent(in) :: mass           ! mass of molecule
        REAL*8             :: ke              ! output
        INTEGER I,J,K
        REAL*8 vel_mid, com_vel
        ke=0.0

        ! Loop over all molecules
        DO I=1,nMol
          ! Loop over all atoms
          DO J=1,nAtomPer
            ! Loop over x,y,z components of velocity
             DO K=1,3
                vel_mid=0.5*(velold(I,J,K)+vel1(I,J,K))
                ke=ke+(vel_mid)**2    ! sum vel^2
            END DO
          END DO
        END DO

        ke = 0.5*ke*mass          ! calculate kinetic energy

      end function

      !-------------------------------------------------------------
      ! function run_SHAKE: "Shake" bonds to fix bond lengths
      !
      ! Input:
      !   natoms:
      !   vel1: (natoms,3) velocities of all atoms
      !   massO: Mass of 1 Ar atom
      !
      ! Output:
      !   ke: Total kinetic energy in the system
      !-------------------------------------------------------------
      function run_SHAKE(nMol,nAPer,pos,pos_old,vel,
     :mass,tbl_SQ,dt) result(mass_term)
        IMPLICIT NONE
        INTEGER, intent(in):: nMol        ! number of molecule
        INTEGER, intent(in):: nAPer       ! number of atoms per mol
        REAL*8:: pos(nMol,nAPer,3)        ! position array
        REAL*8:: pos_old(nMol,nAPer,3)    ! last position array
        REAL*8:: vel(nMol,nAPer,3)        ! velocity array
        REAL*8, intent(in):: mass         ! mass of 1 atom
        REAL*8, intent(in):: tbl_SQ       ! target bond length squared
        REAL*8, intent(in):: dt           ! delta t (time step)

        INTEGER I,J,JA,JB,K,dev_count
        REAL*8 thisLength_SQ,dot_prod,g,correction
        REAL*8 bond_vector(3), old_bond_vector(3)
        REAL*8 mass_term, tau
        REAL*8 tolerance
        LOGICAL flag

        tolerance = 1.0E-7
        mass_term = ((1.0/mass)+(1.0/mass))
        !print*,mass_term

        ! Loop over all molecules
        DO I=1,nMol
          dev_count=0
          flag=.TRUE.

          DO WHILE(flag) ! shake
            flag=.FALSE.

            ! Loop over all atom pairs in molecule
            DO JA=1,nAPer-1
               DO JB=JA+1,nAPer

                  ! Calculate current bond length
                  thisLength_SQ=0.0
                  dot_prod=0.0
                  DO K=1,3
                     bond_vector(K)=pos(I,JA,K)-pos(I,JB,K)
                     thisLength_SQ=thisLength_SQ+bond_vector(K)**2
                     old_bond_vector(K)=pos_old(I,JA,K)-pos_old(I,JB,K)
                     dot_prod=dot_prod+bond_vector(K)*old_bond_vector(K)
                  END DO

                  g=tbl_SQ-thisLength_SQ

                  g=g/(2.0*mass_term*dot_prod)

                  ! Correct bond length (and calculate again)
                  thisLength_SQ=0.0
                  DO K=1,3
                     correction=g*old_bond_vector(K)/mass
                     pos(I,JA,K)=pos(I,JA,K)+correction
                     pos(I,JB,K)=pos(I,JB,K)-correction
                     bond_vector(K)=pos(I,JA,K)-pos(I,JB,K)
                     thisLength_SQ=thisLength_SQ+bond_vector(K)**2
                  END DO

                  ! See if within tolerance
                 IF (thisLength_SQ.NE.0.0) THEN
                    tau=abs(thisLength_SQ-tbl_SQ)/(2.0*tbl_SQ)
                    IF (tau.GE.tolerance) THEN
                       flag=.TRUE.
                    END IF
                    dev_count=dev_count+1
                    IF (dev_count.EQ.500) THEN
                       WRITE(100,*)"too many iterations in SHAKE"
                       WRITE(100,*)'in step ',dev_count
                       WRITE(100,*)'in molecule ',I
                       STOP
                    END IF
                  END IF

                END DO ! JB loop
              END DO ! JA loop

              ! Check if we've met the tolerance
              tau=abs(thisLength_SQ-tbl_SQ)/(2.0*tbl_SQ)
              IF (tau > tolerance) THEN
                 flag=.TRUE.
              END IF
              
              IF (dev_count.EQ.500) THEN
                 WRITE(100,*)"too many iterations in SHAKE"
                 WRITE(100,*)'in step ',dev_count
                 WRITE(100,*)'in molecule ',I
              END IF
              


            END DO ! do while loop

            ! recalculating velocities according to changes in shake
            ! positions

            DO J=1,nAPer
              DO K=1,3
                vel(I,J,K)=(pos(I,J,K)-pos_old(I,J,K))/dt
              END DO
            END DO  !J loop

          END DO ! I loop

        end function


      !-------------------------------------------------------------
      ! Oxgen simulation
      !
      ! Simulate a set of O2 molecules
      ! Based on code from Stephanie Kearing
      !
      ! Input:
      !   Argon_864_initial_coordinates.gro : file with atomic coordinates and
      !       velocities for to initiate simulation
      !
      ! Ouput:
      !   A bunch of files
      !-------------------------------------------------------------


      PROGRAM Oxygen_simulation

      ! Based on code from Stephanie Kearing
      IMPLICIT NONE

      ! Declare functions
      REAL*8 kinetic_energy
      REAL*8 run_SHAKE
      REAL*8 shift_vel_com

      INTEGER I,J,K,L,M,N,time_loop

      ! Set to NVE or NVT
      INTEGER is_NVT ! 0 if no thermostat (i.e. NVE); 1 for NVT

      ! Initiate constants
      REAL*8 NA
      PARAMETER (NA=6.02214076D23)

      REAL*8 BOLTZ
      PARAMETER (BOLTZ=1.38064852D-23) ! in J/K

      ! Set constant parameters for our system
      INTEGER nMol                  ! the number of atoms in the simulation
      PARAMETER (nMol=500)

      INTEGER nAtomPer              ! the number of atoms per molecule
      PARAMETER (nAtomPer=2)

      REAL*8 eps                    ! epsilon value for calculating LJ potential
      PARAMETER (eps=120.0*BOLTZ)   ! epsilon in J

      REAL*8 sigma                  ! sigma value for calculating LJ potential
      PARAMETER (sigma=0.3006000E-9) ! sigma in m
      
      REAL*8 sigma2   
      PARAMETER (sigma2=sigma**2)

      REAL*8 massO                  ! mass of one atom in kg
      PARAMETER (massO=0.016/NA)  

      REAL*8 massO2                 ! mass of one molecule in kg
      PARAMETER (massO2=massO*2.0) 

      INTEGER nstep                 ! number of steps in the simulation
      PARAMETER (nstep=50000)

      INTEGER nsave                 ! frequency to save data
      PARAMETER (nsave=10)

      INTEGER comShiftTime              ! frequency to shift com of velocity
      PARAMETER (comShiftTime=110)

      REAL*8 DT                     ! time step in seconds
      PARAMETER (DT=4.0E-15)

      REAL*8 cutoff,cutoffSQ

      REAL*8 Tref                   ! the reference temprerature in K
      PARAMETER (Tref=320.00000)   

      REAL*8 Length                 ! Length of the box

      REAL*8 bond_length            ! Bond length in m from Perng01
      PARAMETER (bond_length=0.1208E-9) 
      REAL*8 blSQ                   ! Bond length squared
      PARAMETER (blSQ=bond_length**2)

      ! Create arrays for coordinates and velocities
      ! dimensions are: molecules, atoms, (x,y,z)
      REAL*8 pos(nMol,nAtomPer,3),vel(nMol,nAtomPer,3)
      REAL*8 vel_old(nMol,nAtomPer,3),vel_half(nMol,nAtomPer,3)
      REAL*8 pos_old(nMol,nAtomPer,3)

      ! Variables for calculating VLJ and FLJ in x,y,z directions
      REAL*8 rij(3)           ! the vector between atoms I and J
      REAL*8 dist_ij          ! total distance between atoms I and J
      REAL*8 V_ij             ! The LJ potential between atoms I and J
      REAL*8 V_tot            ! The total potential energy in the system
      REAL*8 F_ij             ! the magnitude of force between atoms I and J
      REAL*8 sr_12            ! (sigma/r)^12
      REAL*8 sr_6             ! (sigma/r)^6
      REAL*8 half_box         ! half the length of the box
      
      REAL*8 eps4
      PARAMETER (eps4 = 4.0*eps)
      REAL*8 eps24
      PARAMETER (eps24=24.0*eps)
      
      INTEGER nDOF            ! Atomic DOF
      PARAMETER (nDOF=5)
      REAL*8 KE_Temp          ! Conversion between KE and Temp
      PARAMETER (KE_Temp = 2.0/(REAL(nDOF)*BOLTZ*(nMol)))

      ! Variables for g(r)
      INTEGER nGrBins,ibin
      PARAMETER (nGrBins=500)
      REAL*8 bin_ends(nGrBins)            ! the edge value of each bin for g(r)
      REAL*8 g_of_r(nGrBins)              ! g(r) data
      REAL*8 gr_bin_W,dmax
      PARAMETER (gr_bin_W=0.05E-10,dmax=REAL(nGrBins)*gr_bin_W)

      ! Variables to attempt leap frog generation of new r,v,a
      REAL*8 force(nMol,nAtomPer,3)   ! the force on each particle
      REAL*8 this_force(3)            ! force on this particle
      REAL*8 accel(nMol,nAtomPer,3)   ! acceleration
      REAL*8 is_NVT_scale,KE,vel_2    ! for KE calculation

      REAL*8 time ! the step multiplied by DT to give the time in sec
      REAL*8 temp ! the temperature of the system to be calc from KE

      ! preventing the flying ice cube
      REAL*8 vel_sum(3),vel_scale(3)

      ! reading in/values from input file
      INTEGER idum ! loop/array variables / dummy integer
      CHARACTER*5 resname,atomname1,atomname2

      ! File for storing energy data
      OPEN(92,FILE="Energies_O2_320.txt")
      OPEN(91,FILE="O2_positions_320.gro") ! open output traj file
      OPEN(94,FILE="Temp_320.txt") ! open output temperature
      OPEN(95,FILE="Velocity_O2_320.txt")
      OPEN(96,FILE="Velocity_vector_O2_320.txt")
      OPEN(97,FILE="G_r_all_distances_O2_320.txt")
      OPEN(98,FILE="Oxygen_NVT_320.gro")

      !Set to NVT or NVE
      is_NVT=1

      cutoff=100*sigma
      cutoffSQ=cutoff**2 ! cutoffs for energy equations

      !----------------------
      ! Read in intial file
      !----------------------
      OPEN(11,FILE="Oxygen_500_initial_coordinates.gro") ! open input file
      READ(11,*)
      READ(11,*)
      DO I=1,nMol
          READ(11,31)idum,resname,atomname1,idum,
     :(pos(I,1,K),K=1,3),(vel(I,1,K),K=1,3)
          READ(11,31)idum,resname,atomname2,idum,
     :(pos(I,2,K),K=1,3),(vel(I,2,K),K=1,3)

      END DO
      READ(11,*)Length,Length,Length
      CLOSE(11)

      !--------------------------------------------
      ! unit conversion for pos/vel/box read in from the gro file; to SI
      ! Also calculate average velocity and initialize some variables
      !--------------------------------------------
      vel_sum(1) = 0
      vel_sum(2) = 0
      vel_sum(3) = 0
      Length=Length*1.0E-9 ! from nm to m

      ! Loop over all atoms
      DO I=1,nMol
        DO J=1,nAtomPer
          DO K=1,3                       ! Loop over x,y,z components
              pos(I,J,K)=pos(I,J,K)*1.0E-9    ! from nm to m
              vel(I,J,K)=vel(I,J,K)*1.0E3     ! from nm/ps to m/s
              accel(I,J,K) = 0.0            ! initialize acceleration to 0
              force(I,J,K) = 0.0            ! initialize forces to 0
              pos_old(I,J,K) = pos(I,J,K)
            END DO
         END DO
      END DO
      half_box = Length/2.0D0


      ! Calculate kinetic energy and temperature
      KE = kinetic_energy(nMol, nAtomPer, vel,vel,massO)
      temp = (KE*KE_Temp)
      print*,temp
      is_NVT_scale = sqrt(Tref/temp)  ! Scale for velocity
      print*,is_NVT_scale
      ! Loop over all atoms and scale velocity
      DO I=1,nMol
        DO J=1,nAtomPer
          DO K=1,3
            vel(I,J,K)=vel(I,J,K)*is_NVT_scale
          end do
        end do
      end do

      ! Prevent flying ice cube
      !vel_scale(1) = shift_vel_com(nMol,nAtomPer, vel)

      ! Calculate kinetic energy and temperature
      KE = kinetic_energy(nMol, nAtomPer, vel,vel,massO)
      temp = (KE*KE_Temp)

      DO I=1,nGrBins
         g_of_r(i) = 0.0
         bin_ends(i) = (gr_bin_W*I)
       END DO
       
      !--------------------------------------------
      ! Progress system in time
      !--------------------------------------------
      time = 0.0
      DO time_loop=1,nstep
        time = time+dt                  ! Calculate current time
        V_tot = 0.0                     ! Initialize total velocity to 0

        ! loop over all pairs of atoms
        DO I=1,nMol-1
          DO J=i+1,nMol
            DO M=1,nAtomPer                  ! Atoms in molecule 1
              DO N=1,nAtomPer               ! Atoms in molecule 2

                dist_ij = 0.0             ! reset dist_ij to 0

               
                ! get distance vector between atoms
                DO K=1,3 ! Loop over x,y,z components
                  rij(K)= pos(I,M,K) - pos(J,N,K)
                  rij(K)= rij(K) - length*ANINT(rij(K)/Length)
                  dist_ij = dist_ij + rij(K)*rij(K)
                END DO ! end K: calcualte dist_ij

                ! If within a certain distance, calculate force
                if (dist_ij <= cutoffSQ) then

                  sr_6 = (sigma2/dist_ij)**3
                  sr_12 = sr_6**2

                  ! Calculate potential energy and add to total of system
                  V_ij=eps4*(sr_12-sr_6)
                  V_tot=V_tot+V_ij

                  ! Calculate magnitude of force
                  F_ij=-eps24*(-2.0*sr_12+sr_6)/dist_ij

                  ! Find unit vector of rij and use to add force vector to i and j
                  DO K=1,3
                    this_force(K) = F_ij*rij(K)
                    force(i,M,K) = force(i,M,K) + this_force(K)
                    force(j,N,K) = force(j,N,K) - this_force(K)
                  END DO ! end K: force vectors
                END if  ! end if dist_ij < cutoffSQ
C                 ! Add to g(r) data if in save time
C                 IF(MOD(time_loop,nsave).EQ.0) THEN
C                   dist_ij = sqrt(dist_ij)
C                   ibin=FLOOR((dist_ij)/gr_bin_W)+1
C                   IF (ibin.LE.nGrBins) THEN
C                     g_of_r(ibin)=g_of_r(ibin)+2
C                   END IF
C                 end IF ! if save_loop
              end do ! end N
            end do  ! end M
          END do    ! end J
        END DO      ! end I

        !--------------------------------------------
        ! Use the force to update the positions
        ! This is done using leap frog algorithm
        ! Remember that vel is really vel at t-0.5dt
        !--------------------------------------------

        DO I=1,nMol
          DO J=1,nAtomPer
           DO K=1,3
             pos_old(I,J,K) = pos(I,J,K)
             accel(I,J,K) =force(I,J,K)/massO                  ! Acceleration based on force
             vel_old(I,J,K) = vel(I,J,K)                        ! Store old velocity
             vel(I,J,K)=vel(I,J,K) + DT*accel(I,J,K)              ! Get new velocity (at t+0.5dt)
             pos(I,J,K) = pos(I,J,K) + DT*vel(I,J,K)              ! Get new position (at t+dt)
             vel_half(I,J,K) = (vel(I,J,K) + vel_old(I,J,K))*0.5  ! Get velocity at (t+dt) for KE
             force(I,J,K) = 0.0                               ! reset force array for next time loop
           END DO
          end do

        END DO

        ! Run SHAKE
        ! doesn't use vel, but changes vel and pos
        temp = run_SHAKE(nMol,nAtomPer,pos,pos_old,vel,massO,blSQ,dt)

        ! Calculate kinetic energy and temperature
        KE = kinetic_energy(nMol,nAtomPer,vel_old,vel,massO)
        temp = (KE*KE_Temp)
        WRITE(94,*)temp

        !--------------------------------------------
        ! If we're in NVT, scale velocity to keep set temperature
        !--------------------------------------------
        IF (is_NVT == 1) THEN
            KE = kinetic_energy(nMol,nAtomPer,vel,vel,massO)
            is_NVT_scale = sqrt(Tref/temp)  ! Scale for velocity

            ! Loop over all atoms and scale velocity
            DO I=1,nMol
              DO J=1,nAtomPer
                DO K=1,3
                  vel(I,J,K)=vel(I,J,K)*is_NVT_scale
                end do
              end do
            end do


          !---------------------------------------------------------------------
          ! At set time points, shift com of veloctiy to be 0
          !---------------------------------------------------------------------
          IF(MOD(time_loop,comShiftTime).EQ.0) THEN
            vel_sum(1) = 0.0
            vel_sum(2) = 0.0
            vel_sum(3) = 0.0

            ! Loop over all atoms and find center of mass of velocties
            DO I=1,nMol
              DO J=1,nAtomPer
               DO K=1,3                           ! Loop over x,y,z components
                  vel_sum(K)=vel_sum(K)+vel(I,J,K)  ! Calculate total velocity of system
               END DO
             END DO
            END DO

            ! Calculate how much system needs to shift to keep from flying ice cube
            vel_sum(1) = vel_sum(1)/(nAtomPer*nMol)
            vel_sum(2) = vel_sum(2)/(nAtomPer*nMol)
            vel_sum(3) = vel_sum(3)/(nAtomPer*nMol)
            print*,vel_sum(1)
            ! shift velocity to keep whole system from moving in a direction
            DO I=1,nMol                         ! Loop over all atoms
              DO J=1,nAtomPer
               DO K=1,3                           ! Loop over x,y,z components
                 vel(I,J,K) = vel(I,J,K)-vel_sum(K) ! shift velocity
               END DO
              END DO
            END DO
          END IF
        end if

        !--------------------------------------------
        ! Store system data if we're in a save time
        !--------------------------------------------
        IF(MOD(time_loop,nsave).EQ.0) THEN
           resname='   O'
           IF(is_NVT <2) THEN
             WRITE(91,*)'After step ',time_loop
             WRITE(91,*)(nMol*nAtomPer)

             DO I=1,nMol
                WRITE(91,31)I,resname,atomname1,I,
     :(pos(I,1,K)*1.0E9,K=1,3),
     :(vel(I,1,K)*1.0E-3,K=1,3)
                 WRITE(91,31)I,resname,atomname2,I,
     :(pos(I,2,K)*1.0E9,K=1,3),
     :(vel(I,2,K)*1.0E-3,K=1,3)
              END DO

            WRITE(91,*)Length*1.0E9,Length*1.0E9,Length*1.0E9
            END IF
        END IF
        WRITE(92,*)time,V_tot,KE,V_tot+KE,temp
        
      END DO

      
      WRITE(95,*)'After step ',time_loop
      WRITE(95,*)(nMol*nAtomPer)

      DO I=1,nMol
         WRITE(95,31)I,resname,atomname1,I,
     :(vel(I,1,K)*1.0E-3,K=1,3),
     :(vel(I,1,K)*1.0E-4,K=1,3)
          WRITE(95,31)I,resname,atomname2,I,
     :(vel(I,2,K)*1.0E-3,K=1,3),
     :(vel(I,2,K)*1.0E-4,K=1,3)
       END DO

       WRITE(91,*)Length*1.0E9,Length*1.0E9,Length*1.0E9

      ! Write g_of_r data
      DO I=1,nGrBins
        WRITE(97,*)bin_ends(I),g_of_r(I)
      END DO

      ! Write final coordinates
      WRITE(98,*)'A box of liquid O2'
      WRITE(98,*)(nMol*nAtomPer)
      DO I=1,nMol
        WRITE(98,31)I,resname,atomname1,I,
     :(pos(I,1,K)*1.0E9,K=1,3),
     :(vel(I,1,K)*1.0E-3,K=1,3)
        WRITE(98,31)I,resname,atomname2,I,
     :(pos(I,2,K)*1.0E9,K=1,3),
     :(vel(I,2,K)*1.0E-3,K=1,3)
      END DO
      WRITE(98,*)Length*1.0E9,Length*1.0E9,Length*1.0E9

  31  FORMAT(i5,2a5,i5,3f8.4,3f8.4) ! format of input file
      print*,Length
      CLOSE(92)
      CLOSE(91)
      CLOSE(94)
      CLOSE(95)
      CLOSE(96)
      CLOSE(97)
      END
