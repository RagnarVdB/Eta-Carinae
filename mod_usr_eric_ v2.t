! wind_eta_carinae
module mod_usr
  use mod_hd
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_source        => usrsource
    usr_refine_grid   => usrrefinegrid

    !Code Change to polar
    call set_coordinate_system("polar")
    !call set_coordinate_system("cylindrical")
    !call set_coordinate_system("Cartesian_2D")

    unit_length        = 1.d16 ! cm
    unit_temperature   = 1.d2 ! K
    unit_numberdensity = 2.d0 ! cm^-3

    call hd_activate()

    unit_mass=unit_density*(unit_length**3)

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize variables on the grid
    use mod_physics
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: rad0,xc^D ! rad0: radius of location with wind entry, xc1 xc2 are the coords of the centre
    double precision :: SolarMass=10.d33 ! in grams
    double precision :: PiDouble=3.14159265359 ! Pi
    double precision :: Tbkgrnd, T0A, kT, valtemp(ixI^S)
    double precision :: rhobkgrnd, rho0A, krho, MassPerTime, MassPerYear
    double precision :: velbkgrnd, vel0A, kvel
    logical, save:: first=.true., block0=.true.

    ! check units of variables
    if (mype==0 .and. block0) then
      print *, 'length', unit_length
      print *, 't', unit_time
      print *, 'rho', unit_density
      print *, 'rho_num', unit_numberdensity
      print *, 'mass', unit_mass
      print *, 'p', unit_pressure
      print *, 'T', unit_temperature
      print *, 'v', unit_velocity
      block0=.false.
    endif


    ! outputs message saying experiment type on first program loop iteration
    !xc1=(xprobmin1+xprobmax1)*0.5d0     ! set location of wind

    !Eric's Code
    xc1 = xprobmin1 !The origin should be r=0
    xc2 = xprobmin2 !And theta = 0

    !xc2=(xprobmin2+xprobmax2)*0.5d0     ! set location of wind
    rad0=1.0d16/unit_length             ! inner boundary          =10^16 cm

    T0A=1.0d5/unit_temperature           ! inner source temperature       =10^5 K
    Tbkgrnd=1.d2/unit_temperature        ! outer space temperature        =100  K
    vel0A=2.50d2*1.d5/unit_velocity      ! inner source velocity (part A) =250.0km/s, convert cms, then units
    velbkgrnd=1.00d2*1.d5/unit_velocity ! outer space velocity           =100  km/s, convert cms, then units

    ! rho
    ! rho for correct mass flux
    ! rho = mass flux / (vel * length of border)
    MassPerYear=1.0d-3
    MassPerTime=SolarMass*MassPerYear                                       ! mass loss = 10^-3 solar masses (g) per year
    MassPerTime=MassPerTime*unit_time/(unit_mass*365.25*24*60*60)           ! unit mass per unit time
    rho0A=MassPerTime/(vel0A*2.d0*PiDouble*rad0)                            ! inner source mass density (part A)
    ! Time point for ramping up the mass ejection
    ! rho for outer space is 1 hydrogen masses
    rhobkgrnd=1.d0*1.67d-24/unit_density                                   ! outer space density 10-3 per cm -3


    ! initialise variables with zero value
    w(ixO^S,rho_)=0.0d0
    w(ixO^S,mom(1))=0.0d0
    w(ixO^S,mom(2))=0.0d0
    w(ixO^S,e_)=0.0d0

    !
    ! background density, could use krho/r^2 + rhobckgrnd
    w(ixO^S,rho_)=rhobkgrnd
    ! For no intial jump at boundary: krho=(rho0-rhobkgrnd)*rad0**2
    krho=(rho0A-rhobkgrnd)*(rad0**2)
    !Eric's Change
    !w(ixO^S,rho_)=( krho/({^D&(x(ixO^S,^D)-xc^D)**2+}) + rhobkgrnd)
    !Turned following line into code since divide by 0
    !w(ixO^S,rho_) = (krho/(x(ixO^S,1)-xc1)**2 +rhobkgrnd) !Note that this decays as 1/r^2

    !
    ! background velocity, could use an initial setup like kvel/r^2 + velbkgrnd
    ! this is a little more complicated in Cartesians because we need the correct component for each point in space
    w(ixO^S,mom(:))=0.d0
    kvel=(vel0A-velbkgrnd)*(rad0**2)
    ! background component
    !Eric's Change
    !w(ixO^S,mom(1))= velbkgrnd*(x(ixO^S,1)-xc1)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !w(ixO^S,mom(2))= velbkgrnd*(x(ixO^S,2)-xc2)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    w(ixO^S,mom(1)) = velbkgrnd !the radial velocity is the background velocity
    w(ixO^S, mom(2)) = 0.d0 !the angular velocity is zero
    ! decaying component from central value
    !w(ixO^S,mom(1))=w(ixO^S,mom(1)) + ( kvel/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+}))*(x(ixO^S,1)-xc1)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !w(ixO^S,mom(2))=w(ixO^S,mom(2)) + ( kvel/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+}))*(x(ixO^S,2)-xc2)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !Again, turn following code into comment as to not divide by 0
    !w(ixO^S,mom(1))=w(ixO^S,mom(1)) + kvel/(x(ixO^S,1)-xc1)**2 !Add decay term


    ! background temperature
    ! T =  kT/r^2 + Tbckgrnd
    !Eric's change
    !kT=(T0A-Tbkgrnd)*(rad0)
    kT=(T0A-Tbkgrnd)*(rad0**2) !**2 since planning to decay T as 1/r^2 instead of 1/r
    !w(ixO^S,e_)=( kT/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+}) + Tbkgrnd)
    !Made the following code into comment since if x=xc1 i.e. the origin; the code shall divide by 0
    !w(ixO^S,e_) = (kT/(x(ixO^S,1)-xc1)**2 +Tbkgrnd)

    ! convert T to thermal pressure
    !     T=p/rho, so p=T*rho
    w(ixO^S,e_)= w(ixO^S,e_)*w(ixO^S,rho_)

    ! inside radius
    ! * velocity zero at centre, increasing outward
    ! * set mass uniform inside initially
    !     such that this will produce a rate of mass flux rate equal to the desired amount
    ! * set temperature = T0A = 10^5K via pressure, p=T*rho
    !Eric's Change
    !where(({^D&(x(ixO^S,^D)-xc^D)**2+})<=rad0**2)
    !    w(ixO^S,mom(1))=vel0A*({^D&(x(ixO^S,^D)-xc^D)**2+})/(rad0**2)*(x(ixO^S,1)-xc1)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !    w(ixO^S,mom(2))=vel0A*({^D&(x(ixO^S,^D)-xc^D)**2+})/(rad0**2)*(x(ixO^S,2)-xc2)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !    w(ixO^S,rho_)=rho0A
    !    w(ixO^S,e_)=T0A*w(ixO^S,rho_)
    !endwhere  
    where(((x(ixO^S,1)-xc1))<=rad0)
        !In v2 try to make inside of radius0, a square increase of velocity instead of linear
        w(ixO^S,mom(1))=vel0A*(x(ixO^S,1)-xc1)**2/(rad0)**2
        w(ixO^S,rho_)=rho0A
        w(ixO^S,e_)=T0A*w(ixO^S,rho_)
    endwhere
    !Change from v1: made it so the code doesn't divide by zero
    where(((x(ixO^S,1)-xc1))>rad0)
        !Density
        w(ixO^S,rho_) = (krho/(x(ixO^S,1)-xc1)**2 +rhobkgrnd) !Note that this decays as 1/r^2
        !Velocity
        w(ixO^S,mom(1))=w(ixO^S,mom(1)) + kvel/(x(ixO^S,1)-xc1)**2 !Add decay term
        !Temperature
        w(ixO^S,e_) = (kT/(x(ixO^S,1)-xc1)**2 +Tbkgrnd)
        !Convert into Thermal pressure
        w(ixO^S,e_)=T0A*w(ixO^S,rho_)
    endwhere  


    if (first) then
       if (mype==0) then
          print *, 'Source radius, rad0', rad0
          print *, 'Background temperature', Tbkgrnd
          print *, 'Part A source temperature', T0A
          print *, 'Part A source velocity', Vel0A
          print *, 'Part A source density', rho0A
          print *, 'Part A solar masses per year', MassPerYear
          print *,'2D HD wind in Cartersian coordinates'
       end if
       first=.false.
    end if


    ! as set pressure, velocity, convert these to energy and momentum (the code's conserved values)
    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_physics
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: tmp(ixI^S) 

    ! You can add output to the "converted" files (e.g. the non .dat ones)
    ! Temperature can be calculated in post processing via p/rho
    ! Here it is done automatically by adding it to the list of w variables
    call hd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! output the temperature p/rho
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te'

  end subroutine specialvarnames_output

  subroutine usrsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
  ! Purpose:
  !   Add in wind terms from the eta carinae system at the inner boundary
  !   Called each step
  !     Gonzalez 2022
  !     originates at a radius of r_0 = 10^16 cm
  !
  ! Parameter list
  !   qdt :  heating "timestep" delta t_Q
  !   ixI^l: limits (min max) of inner cell indices (no ghost cells)
  !   ixO^l: limits (min max) of outer cell indices (inc ghost cells)
  !   iw^l:  limits of variable indices
  !   ?qtC : some energy terms? handed to source subroutines
  !   ?wCT : some variable terms handed to source subroutines
  !   qt :   experiment time
  !   w:     variable values
  !   x:     Cell coords
    use mod_physics

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rad0,xc^D ! rad0: radius of location with wind entry, xc1 xc2 are the coords of the centre
    double precision :: SolarMass=10.d33 ! in grams
    double precision :: PiDouble=3.14159265359 ! Pi
    double precision :: T0A, valtemp(ixI^S)
    double precision :: rhobkgrnd, rho0A, MassPerTime
    double precision :: velbkgrnd, vel0A, temp
    double precision :: timepoint1
    logical, save:: first=.true.

    !Eric's change
    !xc1=(xprobmin1+xprobmax1)*0.5d0      ! set location of wind
    !xc2=(xprobmin2+xprobmax2)*0.5d0      ! set location of wind
    xc1 = xprobmin1 !The origin should be r=0
    xc2 = xprobmin2 !And theta = 0

    rad0=1.0d16/unit_length              ! inner boundary          =10^16 cm
    T0A=1.0d5/unit_temperature           ! inner temperature       =10^5 K
    vel0A=2.50d2*1.d5/unit_velocity      ! inner velocity A        =250   km/s, convert cms, then units

    ! rho
    !
    ! rho for correct mass flux
    ! rho = mass flux / (vel * length of border)
    MassPerTime=SolarMass*1.0d-3                                            ! mass loss = 10^-3 solar masses (g) per year
    MassPerTime=MassPerTime*unit_time/(unit_mass*365.25*24*60*60)           ! unit mass per unit time
    rho0A=MassPerTime/(vel0A*2.d0*PiDouble*rad0)                            ! inner mass density A
    rhobkgrnd=1.d0*1.67d-24/unit_density                                   ! outer space density 10-3 per cm -3
    !velbkgrnd=0.0d0


    ! This was an old version where I used timepoint parameter to increase some values up slowly, decided against it!
    !timepoint1=1.0d1*(365.25*24*60*60)/unit_time                             ! ramp up the rho over 10 years
    !if(qt .lt. timepoint1) then
    !    rho0A=rhobkgrnd+(rho0A-rhobkgrnd)*(qt/timepoint1)
    !endif
    !if(qt .lt. timepoint1) then
    !    vel0A=velbkgrnd+(vel0A-velbkgrnd)*(qt/timepoint1)
    !endif


    ! convert into velocity etc, not needed as setting everything here
    call hd_to_primitive(ixI^L,ixO^L,w,x)

    ! inside radius
    ! * velocity zero at centre, increasing outward
    ! * set mass uniform inside
    !     such that this will produce a rate of mass flux rate equal to the desired amount
    ! * set temperature = T0A = 10^5K via pressure, p=T*rho
    !Eric's Change
    !where(({^D&(x(ixO^S,^D)-xc^D)**2+})<=rad0**2)
    !    w(ixO^S,mom(1))=vel0A*({^D&(x(ixO^S,^D)-xc^D)**2+})/(rad0**2)*(x(ixO^S,1)-xc1)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !    w(ixO^S,mom(2))=vel0A*({^D&(x(ixO^S,^D)-xc^D)**2+})/(rad0**2)*(x(ixO^S,2)-xc2)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !    w(ixO^S,rho_)=rho0A
    !    w(ixO^S,e_)=T0A*w(ixO^S,rho_)
    !endwhere  
    where(((x(ixO^S,1)-xc1))<=rad0)
        !Change from v1: square velocity instead of linear
        w(ixO^S,mom(1))=vel0A*(x(ixO^S,1)-xc1)**2/(rad0)**2
        w(ixO^S,rho_)=rho0A
        w(ixO^S,e_)=T0A*w(ixO^S,rho_)
    endwhere  

    ! as set pressure, velocity, convert these to energy and momentum (the code's conserved values)
    call hd_to_conserved(ixI^L,ixO^L,w,x)

    !Add Extra wave from cartesian
    !if (qt > 1.0d-1) then
    !    MassPerTime=SolarMass*0.7d0                                           
    !    MassPerTime=MassPerTime*unit_time/(unit_mass*365.25*24*60*60)
    !    vel0A=1.0d2*3.d5/unit_velocity      ! inner velocity A        =250   km/s, convert cms, then units

    !    rho0A=MassPerTime/(vel0A*2.d0*PiDouble*rad0)                            ! inner mass density A
    !    rhobkgrnd=1.d0*1.67d-24/unit_density  
    !    T0A=1.0d5/unit_temperature           ! inner temperature       =10^5 K
        
    !    call hd_to_primitive(ixI^L,ixO^L,w,x) 
        
    !    call hd_to_conserved(ixI^L,ixO^L,w,x)
    !    where(({^D&(x(ixO^S,^D)-xc^D)**2+})<=rad0**2)
    !        w(ixO^S,mom(1))=vel0A*({^D&(x(ixO^S,^D)-xc^D)**2+})/(rad0**2)*(x(ixO^S,1)-xc1)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !        w(ixO^S,mom(2))=vel0A*({^D&(x(ixO^S,^D)-xc^D)**2+})/(rad0**2)*(x(ixO^S,2)-xc2)/sqrt({^D&(x(ixO^S,^D)-xc^D)**2+})
    !        w(ixO^S,rho_)=rho0A
    !        w(ixO^S,e_)=T0A*w(ixO^S,rho_)
    !    endwhere  
    !endif
    if (qt > 1.0d-1) then
        MassPerTime=SolarMass*0.7d0                                           
        MassPerTime=MassPerTime*unit_time/(unit_mass*365.25*24*60*60)
        vel0A=1.0d2*3.d5/unit_velocity      ! inner velocity A        =250   km/s, convert cms, then units

        rho0A=MassPerTime/(vel0A*2.d0*PiDouble*rad0)                            ! inner mass density A
        rhobkgrnd=1.d0*1.67d-24/unit_density  
        T0A=1.0d5/unit_temperature           ! inner temperature       =10^5 K
        
        !I don't understand why the call hd primitive and hd conserved are done here, so i'll comment them
        !call hd_to_primitive(ixI^L,ixO^L,w,x) 
        !call hd_to_conserved(ixI^L,ixO^L,w,x)

        !Here converted to polar
        where((x(ixO^S,1)-xc1)<=rad0)
            !Once again I'll try out a square inside of the radius instead of a linear increase
            w(ixO^S,mom(1))=vel0A*(x(ixO^S,1)-xc1)**2/(rad0)**2
            w(ixO^S,rho_)=rho0A
            w(ixO^S,e_)=T0A*w(ixO^S,rho_)
        endwhere  
    endif

    !Added call conserved later; just like Malcolm did 
    call hd_to_conserved(ixI^L,ixO^L,w,x)


  end subroutine usrsource 

  subroutine usrrefinegrid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
    ! Purpose:
    !   Enforce additional refinement or coarsening
    !     ignore iprob 1 option, usr int switch for setting up the model with =1
    !
    ! Parameter list
    !   igrid: main program gives -grid number and current refinement level to inspect
    !   level: main program gives grid number and -current refinement level to inspect
    !   ixi^l: Limits (min max) of inner cell indices (no ghost cells)
    !   ixo^l: Limits (min max) of outer cell indices (inc ghost cells)
    !   qt:    Experiment time
    !   w:     Variable values
    !   x:     Cell coords
    !
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen
    logical :: rfgrid
    double precision :: rad0,xc^D ! rad0: radius of location with wind entry, xc1 xc2 are the coords of the centre

    ! Enforce grid refinement at start of experiment to be hyperfine near the central source
    !Eric's Changes
    !xc1=(xprobmin1+xprobmax1)*0.5d0
    !xc2=(xprobmin2+xprobmax2)*0.5d0
    xc1 = xprobmin1
    xc2 = xprobmin2
    rad0=1.0d16/unit_length
    ! fix the inner layer to refinement level 6
    !Eric's Changes
    !if(qt .le. 1.0d-12) then
    !  if (any({^D&(x(ixO^S,^D)-xc^D)**2+}  <= ((2.0*rad0)**2))) then
    !  if (any({^D&(x(ixO^S,^D)-xc^D)**2+}  >= ((0.5*rad0)**2))) then
    !    if (level<5) then
    !      refine=1
    !      coarsen=-1
    !    else
    !      refine=-1
    !      coarsen=-1
    !    endif
    !  endif
    !endif
    !endif

    if(qt .le. 1.0d-12) then
      if (any((x(ixO^S,1)-xc1)  <= ((2.0*rad0)))) then
      if (any((x(ixO^S,1)-xc1) >= ((0.5*rad0)))) then
        if (level<5) then
          refine=1
          coarsen=-1
        else
          refine=-1
          coarsen=-1
        endif
      endif
    endif
    endif
  end subroutine usrrefinegrid

  ! Just a function for if you are having NAN problems and need to check whether w contains them!
  elemental logical function isNANpresent(var)
    double precision, intent(in) :: var
    isNANpresent=(ieee_is_nan(var))
  end function isNANpresent

end module mod_usr
