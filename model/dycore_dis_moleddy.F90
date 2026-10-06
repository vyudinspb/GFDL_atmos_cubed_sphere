!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief Molecular and eddy (LES) physics-based dissipative processes
!> a)   upper layer molecular viscosity, heat conductivity and
!>     diffusion of minor species (major species cannot be diffused !!! as in WAM-IPE) 
!>
!>       viscosity (mur), conductivity (lam), and diffusivity (d12)
!> 
!>@author H.-M. H. Juang, NOAA/NWS/NCEP/EMC in (May 2020), as original explicit diifusion solver
!>        on the vertically Lagragian layers
!>@Kevin Viner (NOAA/EMC) and Valery Yudin (NOAA/NESDIS/SWO-CUA) (2025) applications
!          on the Eulerian presssure layers rho*Cp*dT/dt.... with sub-cycling
!>@Valery Yudin (NOAA/NESDIS/SWO-CUA) (2026) , adding eddy LES dissipation and restoring
!          energy conservation due to the viscous loss of kinetic energy
!>@Valery Yudin (NOAA/NESDIS/SWO-CUA) (Jan 2026) adding contribution of He in the molecular
!         viscosity and heat conductivity for the top-lid (850 km) UFSWAM configurations
!>@Valery Yudin (NOAA/NESDIS/SWO-CUA) (May 2026) new dycore_dis_moleddy.F90
!                                     DAF-LES/MED configuration of eddy-mol dis-n as a sub-step
!                                     of UFSWAM-dynamical core
!
module dycore_dis_moleddy_mod

! Modules Included:
! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>rdgas, cp_air</td>
!   </tr>
! </table>

      use constants_mod,      only: rdgas, cp_air
      use fv_mp_mod,          only: is_master, mp_reduce_max
      use mpp_mod,            only: FATAL, mpp_error, stdlog, input_nml_file
      use fms_mod,            only: check_nml_error, open_namelist_file, close_file
      use fv_grid_utils_mod,  only: g_sum, cubed_to_latlon, update_dwinds_phys
      use mpp_domains_mod,    only: domain2d
      use fv_arrays_mod,      only: fv_grid_type, fv_grid_bounds_type, fv_flags_type
      use fv_mp_mod,          only: is_master
      use fv_mp_mod,          only: start_group_halo_update, complete_group_halo_update
      use fv_mp_mod,          only: group_halo_update_type
      use mpp_domains_mod,    only: CGRID_NE, DGRID_NE, mpp_get_boundary, mpp_update_domains
      use fv_timing_mod,      only: timing_on, timing_off
      use tp_core_mod,        only: copy_corners, &
                              deln_flux_explm, deln_flux_explm_udvd

      use fv_diagnostics_mod,  only: prt_mxm

#ifdef MULTI_GASES
      use multi_gases_mod,   only: ind_gas, num_gas, virqd, vicpqd, vicvq, vicvqd, virq
#endif


      implicit none
      integer :: ind_gas_str, ind_gas_end      
      
! namelist molecular_diffusion_nml

      integer :: md_layers, md_tadj_layers
      real    :: tau_visc, tau_cond, tau_diff
      real    :: md_consv_te
      logical :: md_subcycle 
!      
! 2026 md_helium flag md_helium=.true.  high-top UFSWAM configuration
!              
      logical :: md_helium
!
!  
! namelist /eddy_les_nml/    
! 
      real :: les_kzmin, les_kzmax, les_khmin, les_khmax ! max-min eddy viscosity of SLES
      real :: les_lzmin, les_lhmin, les_lzmax, les_lhmax ! max-min values of eddy lenghts
      real :: les_zgrowth, les_hgrowth  ! vertical growth par-rs for les-scales with height
      real :: les_prnum                 ! Les_Pr = Kv/Kt 
      real :: les_schnum                ! Les_Schmidt = Kv/Diff
      real :: lturbmax                  ! 1/4*delz
      
      real, parameter :: lturb =33.     ! turbulent scale ~ 1./30 of 1-km 
      
      real, parameter :: lturbmaxC =10.*lturb     ! turbulent scale ~ 1./3 of 1-km  
              
      real, parameter :: dw2min = 1.e-2 ! min(uz*uz+vz*vz)  
      real, parameter :: ricrit =0.25   ! Richardson cr. number
      real, parameter :: bn2min = 1.e-6 ! BV2-min   
          
      logical ::  les_dyn      ! .true./.false. Dynamic LES of Germano et al. & Lilly (under dev-nt)
      integer ::  les_sample   ! number of grid-boxes used in sampling of unresolved scales
      real    ::  les_lheff
      real    ::  les_lzeff  
         
      real, parameter:: Runiv = 8314.4598 
      real, parameter:: amo=15.9994, amo2=31.9999, amo3=47.9982     !g/mol
      real, parameter::              amn2=28.013,  amh2o=18.0154    !g/mol
!==========================================================================
!      
! define some constants A*T^0.69 ~> T^0.69 -dependence og 
!                       of  Banks & Kockarts, 1973 can provide large errors
!   Coefficients  based on Pavlov (2017, 2019) will be tested next
!  
!==========================================================================  
!< muo3 and muh2o are not precise, correct later
! corrected https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=907539
!
      real, parameter:: muo=3.9e-7, muo2=4.03e-7,  muo3=1.49e-7     !kg/m/s
      real, parameter::             mun2=3.43e-7,  muh2o=9.70e-7    !kg/m/s
      
!< lao3 is not precise values, but o3_n is very small

      real, parameter:: lao=75.9e-5, lao2=56.e-5,  lao3=36.e-5     !kg/m/s
      real, parameter::              lan2=56.e-5,  lah2o=55.e-5    !kg/m/s
      real, parameter:: cpo=1299.185, cpo2=918.0969, cpo3=820.2391
      real, parameter::               cpn2=1031.108, cph2o=1846.00
      real, parameter:: avgd=6.0221415e23  ! Avogadro constant
      real, parameter:: bz=1.3806505e-23   ! Boltzmann constant J/K
      real, parameter:: a12=9.69e18 ! O-O2 diffusion params
      real, parameter:: s12=0.774
!
! 2026: V. Yudin: helium He and H for the top-lid Starlink density predictions
!
      real, parameter:: cphe = 5193.1631, muhe = 3.84e-7, lahe= 299.e-5 
      real, parameter:: cphl = 14307.*2,  muh1 =1.22e-7,  lah1 =379.e-5     
      real, parameter:: R_he = 2076.9, Rh1 = 2.*4124.0      
      real, parameter:: amh1 = 1.00784, amhe = 4.002602
!
      public :: edmol_diss_init
      public :: read_namelist_edmol_nml
      public :: eddy_molec_coefs, edcoef_sles
      public :: dycore_disedmol_step     

      CONTAINS
! --------------------------------------------------------
      subroutine edmol_diss_init(ncnst,nwat)
!--------------------------------------------
! edmol_diss_init for each effect from namelist
! Input:
!        ncnst    : number of all prognostic tracers
!        nwat     : number of water and the end location of water
!--------------------------------------------
      integer, intent(in):: ncnst, nwat
!

      if( is_master() ) then
        write(*,*) ' eddy_molec_dissipation is on'
        write(*,*) ' molecular_diffusion number of layers ',md_layers
        write(*,*) ' molecular_diffusion subcycling ',      md_subcycle	
        write(*,*) ' molecular_diffusion for Helium ',      md_helium
      endif

#ifdef MULTI_GASES
      ind_gas_str = ind_gas
      ind_gas_end = num_gas
#else
      ind_gas_str = nwat + 1
      ind_gas_end = ncnst
#endif

      return
      end subroutine edmol_diss_init

! --------------------------------------------------------
      subroutine read_namelist_edmol_nml(nml_filename,ncnst,nwat)

      character(*), intent(IN) :: nml_filename
      integer, intent(IN) :: ncnst, nwat
      integer :: ierr, f_unit, unit, ios

      namelist /molecular_diffusion_nml/ tau_visc, tau_cond, tau_diff, &
                                         md_layers, md_tadj_layers, &
                                         md_consv_te, md_subcycle, md_helium 
					 
      namelist /eddy_les_nml/ les_kzmin, les_kzmax, les_khmin, les_khmax, &
                              les_lzmin, les_lhmin, les_lzmax, les_lhmax, &
			      les_zgrowth, les_hgrowth, les_prnum, les_schnum, les_dyn, les_sample

      tau_visc  = 1         !  not active
      tau_cond  = tau_visc  !  not active
      tau_diff  = tau_visc  !  not active
      md_tadj_layers = 0    !  not active
      md_consv_te =.true.   !  not active
      
      md_helium =.false.
      md_layers = 70        ! active only in old version domain in layers where Hor_MDiF is applied    
      
      unit = stdlog()

! Read molecular_diffusion namelist
      read (input_nml_file,molecular_diffusion_nml,iostat=ios)
      ierr = check_nml_error(ios,'molecular_diffusion_nml')

      write(unit, nml=molecular_diffusion_nml)
      call edmol_diss_init(ncnst,nwat)
!      
!old-spec in PJ_VSP/fv_mapz.F90: real, parameter :: lz_smag = 35., cs_smag = 0.53  
! 
      les_khmax = 1.e7           !m2/s
      les_khmin = 1.             !m2/s   
      les_kzmax = 1.e2           !m2/s
      les_kzmin = 1.e-2          !m2/s  
      les_lzmin = 30.            !meters
      les_zgrowth = 21.e3        ! ~3H in exp(Z/3H)
      les_hgrowth = 200.e3
      les_lzmax = les_lzmin*100.     
      les_lhmin = 1.e3           ! ~1 km
      les_lhmax = les_lhmin*100. ! 100 km ~ C96 resolution       
      les_prnum = 1. 
      les_schnum =1.
      les_dyn =.false.           ! flag for Germano-Lilly DLES
      les_sample = 4             ! 4dx x 4dy  
!
! efficiencies of unresolved eddies in grid-box 1./4 of grid-size
!      
      les_lheff = 1./les_sample/les_sample 
      les_lzeff = 1./les_sample 
      
      read (input_nml_file,eddy_les_nml,iostat=ios)
      ierr = check_nml_error(ios,'eddy_les_nml')
      
      return
      end subroutine read_namelist_edmol_nml

! --------------------------------------------------------
      subroutine eddy_molec_coefs(dim, ispo, ispo2, isphe, plyr, temp, q, mur, lam, d12, cp2, molw)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Input: dim : length of field
!        plyr: pressure at given layer (pascal)
!        temp: temperature (K)
!        q   : tracers, needs only gas (kg/kg)
! Ouput: mur : viscosity for momemtum
!        lam : conductivity for thermodynamics
!        d12 : diffusivity for tracers
!--------------------------------------------
      integer, intent(in):: dim, ispo, ispo2, isphe
      real, intent(in):: plyr(dim), temp(dim), q(dim,*)
      real, intent(out):: mur(dim), lam(dim), d12(dim), cp2(dim), molw(dim)
! Local:
      integer i, n
      real   am, a12bz, avgdbz
      real   qo,  qo2,  qo3,  qn2,  qh2o, qhe
      real   o_n, o2_n, o3_n, n2_n, h2o_n, he_n 
      real   mu, la, t69, rho, cpx, cvx, rho1

!constants
      a12bz = a12  * bz
      avgdbz= avgd * bz

      	
      do n=1,dim
!check
!        if(plyr(n).le.0.0) call mpp_error(FATAL,"ERROR non positive value of plyr")
!        if(temp(n).le.0.0) call mpp_error(FATAL,"ERROR non positive value of temp")

        d12(n) = a12bz*temp(n)**s12 * temp(n)/plyr(n)
	
!        if( d12(n).lt.0.0 .and. is_master() ) then
!          write(*,*) 'ERROR negative d12 a12bz temp plyr s12 ', &
!                      d12(n),a12bz,temp(n),plyr(n),s12
!        endif

        qo   = max(0.0,q(n,ispo ))
        qo2  = max(0.0,q(n,ispo2))
 
	
        if (md_helium) then
	  qhe  = max(0.0,q(n,isphe))	
	endif
	
        qn2 = 1.0 - qo - qo2	
	if (md_helium) qn2 = qn2 -qhe
	
        qn2  = max(min(qn2,1.0),0.0)	
!
! 3 major gases only
!		
        cpx  = cpo*qo + cpo2*qo2 + cpn2*qn2 
	if (md_helium) cpx  =	cpx + qhe * cphe	
	
	cp2(n) = cpx
	 
        am = qo/amo + qo2/amo2 + qn2/amn2
	if (md_helium)  am = qhe/amhe +am
	
        am = 1.0 / am                       ! g/mol [29-16]
	molw(n) =  am
	
	cvx  = cpx - runiv/am
	
        o_n   = qo   * am / amo
        o2_n  = qo2  * am / amo2
        n2_n  = qn2  * am / amn2
	if (md_helium) he_n  = qhe  * am / amhe
	
        mu = o_n*muo + o2_n*muo2 + n2_n*mun2
	if (md_helium) mu = mu +he_n*muhe	
	
        la = o_n*lao + o2_n*lao2 + n2_n*lan2
	if (md_helium) la = la +he_n*lahe
	
        t69 = temp(n) ** 1.69                          ! T^(0.69+1) 1/rho=(R/molw*T)/P
!        rho = 1e-3 * am * plyr(n)/temp(n) / avgdbz    ! kg/m^3   avgd * bz/molw = Runiv/molw =8314./29=287
        rho1 =avgdbz/plyr(n)/am
!	
! kinematic molecular visc/cond Banks and Kockarts, Aeronomy,[1973]
! will be updated following expressions  by Pavlov (2017, 2018) 
!	
        mur(n) = mu * t69 * rho1
!vay-bug        lam(n) = la * t69 / rho / cvx only for control-volume Z-coordinate system
        lam(n) = la * t69 *rho1 / cpx

      enddo
     
      return
      
      end subroutine eddy_molec_coefs
      
  subroutine edcoef_sles(ng, npx, npy, bd, u,  v, w, delz2, vumol, edvis, edthc, edmix, gridstruct) 
!
! here w-cell centered (u, v) on D-grid
!  
  use sw_core_mod,      only :	sles_tens_shear      
  type(fv_grid_bounds_type), intent(IN) :: bd     
  integer, intent(in):: npx, npy, ng
  
  real, intent(inout):: w(bd%isd:bd%ied,bd%jsd:bd%jed)
  real, intent(inout):: u(bd%isd:bd%ied,bd%jsd:bd%jed+1)
  real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed)      
  real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed)  :: vumol
  real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed)  :: delz2    
  type(fv_grid_type),  intent(INOUT), target :: gridstruct 
  
  real,  dimension(bd%isd:bd%ied,bd%jsd:bd%jed), intent(out) ::  edvis, edthc, edmix  
!  
! local
!      real :: les_kzmin, les_kzmax, les_khmin, les_khmax ! max-min eddy viscosity of SLES
!      real :: les_lzmin, les_lhmin, les_lzmax, les_lhmax ! max-min values of eddy lenghts
!      real :: les_zgrowth, les_hgrowth  ! vertical growth par-rs for les-scales with height


   real :: tens_uv(bd%isd:bd%ied,bd%jsd:bd%jed)     !< Smag-tensor strain/shear on sphere
   real :: Lhdis2                                   ! [Lh-length]^2 or Area of the mesh  
   real :: ed_kh, p23         
   integer :: i, j
   integer :: is,  ie,  js,  je  
      p23= 2./3.
      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je   
      
     call sles_tens_shear(u, v, tens_uv, bd, npx, npy, gridstruct, ng) 
!     
!  Approximate the horiz-isotropic Kh-eddy coefficients: edvis, edthc, edmix
!
! vertical wind shear Sw =sqrt[(wx^2+wy^2)] can be added for the NH dycores
! tens_b is defined at         [do j=js,je+1 ; do i=is,ie+1 ]
! Option 1: Lh ~ sqrt(da_min_c or da_min) * [efficiency = 1/4 ]
       Lhdis2 = gridstruct%da_min * les_lheff    ! eddy length x les_Lh_eff ~ 1/(4)
!
! Option 2: Lh modulated by delz-variable "quasi-isotropic" mixing length" 
!           Lhdis2 = (dx * dy * dz)**(2./3) analog of LES
!           on the sphere Lhdis2 = [r*r*dr*cos(lat)*dlon*dlat= A_xy*dr)**(2./3)
!       
       do j=js,je+1
          do i=is,ie+1
	     Lhdis2=(gridstruct%da_min*delz2(i,j))**p23 ! delz-variable L2_mixlen
	     
	     ed_kh = tens_uv(i,j)*Lhdis2                  ! tensorH * Lh2
             ed_kh = min(ed_kh, les_khmax)                ! maxval for hor-eddy mixing
	     
	     edvis(i,j) = max(0.01*vumol(i,j), ed_kh) ! if ed_kh < 0.01*vumol => add 10% of vumol 
                                                      ! min(ed_kh) = 0.01*vumol(i,j) in the thermosphere	     
          enddo
       enddo  
        
       edthc = edvis/les_prnum         ! les_prnum = les_schnum = 1 or  3 as for UGWP-v1 physics
       edmix = edvis/les_schnum        !
                      
   end subroutine edcoef_sles   
!   
   subroutine get_linear_diss(npz, is, ie, t2, z2, kion, radnc)
   
   implicit none
!
! temporary specification of the linear dissipation 
!   newtonian cooling =const ;  ion drag height-dep in the ionosphere
!   
   integer, intent(in) :: npz, is, ie
   real, intent(in)      ::   t2(is:ie, npz), z2(is:ie, npz)
   real, intent(out)     :: kion(is:ie, npz), radnc(is:ie, npz)
   
!local   
   integer :: k, j, i ,ii, kk
   integer :: ktrop  
   real :: Ne_prof, nef2, hef2, hch_f2, nede, hkde, hpde, bdde, bede,zchap
   
    ktrop = 90                      ! ktrop the tropopause level from the surface
    do k= 1, npz
     do i=is, ie
     radnc(i, k) = 1./365./86400.   ! newtonian cooling inverse time scale (1/sec) now min-value 1-Year
                                    ! global-tabulated values aka WACCM for GW-damping
     kion(i, k) =  0.               ! zero ion-drag proxi analytical expression requires e-dens profile see below
     enddo
    enddo
    
   return
!
! Simple Ne_prof (Chapman-layer) & kion =function (T, Ne, Z)
! 
    nef2 = 7.10264e5
    hef2 = 300.e3
    hch_f2 = 39.e3
    nede = 1.43e13
!    
    hkde = 95.e3
    hpde = 70.e3
    bdde = 0.65e-3
    bede =0.14e-3 
!    
   do ii = is, ie
   do kk = npz, ktrop,-1  
      zchap = (z2(ii,kk)-hef2)/hch_f2 
      Ne_prof = nef2*exp(0.5*(1- zchap -exp(-zchap)))
      kion(ii,kk) = 7.22e-11* (t2(ii, kk) )** 0.37 *Ne_prof 
   enddo 
   enddo	 
   end subroutine get_linear_diss
!   
   subroutine edcoef_sles_vert(jind, je, is, ie, npz, t2, u2, v2, dp2, dz2, zm2, ht2,    &
                               pmid, pint, pmid_u, pint_u, pmid_v, pint_v, cp2, grav2,   &
	                       kap2, vu2, kt2, df2, kzw, kzu, kzt, kzq, kzv, con_adj)
!
! defines efffective Kz-coefficents of eddy viscosity, conductivity and difuusion of minor species
!         major species that define air density are not diffused (He-O-O2-N2) 
!
   integer, parameter  :: flag_lzmix = 1  ! (1): for delz/4* exp(Zgrow); (2) for analyt-profile lzmix	
   real,    parameter  :: lz_surf = 35. , lz_top = 110., lz_turb = 75.,  zm_turb =100.e3    
   integer, intent(in) :: jind, je
   integer, intent(in) :: is, ie, npz
   real, intent(in), dimension(is:ie, npz)   :: v2, u2, dp2, dz2, pmid, pmid_u, cp2, vu2, kt2, df2
   real, intent(in), dimension(is:ie, npz)   :: zm2, ht2, kap2   
   real, intent(in), dimension(is:ie, npz+1) :: pint, pint_u, grav2
   real, intent(in), dimension(is:ie+1, npz+1) :: pint_v
   real, intent(in), dimension(is:ie+1, npz)   :: pmid_v
   real, intent(inout), dimension(is:ie, npz)   :: t2      
   real, intent(out), dimension(is:ie, npz+1)    :: kzt, kzq 
   real, intent(out), dimension(is:ie+1, npz+1)  :: kzv
   real, intent(inout),dimension(is:ie, npz+1)   :: kzu, kzw  
   logical :: con_adj 
!
! local
!   	
   integer :: i, j, k, ik, ii, kk
   real  :: ht_up, ht_bot, ztop, zbot
   real  :: grav(npz+1),    pe(npz+1), kzz_edd(npz+1)     
   real  :: t1(npz), pt(npz), tadj(npz)
   real  :: zgrow, gmid
   real  :: zlturb     
   real  :: uz, vz, tz, tc, bn2, dzg, gcp, rik, shr2, ritur
   real  :: bn_uns, zmetk, w1, kamp, dtemp, ht    
   integer  :: ins_con(npz), num_ins 
   integer  :: nfixed, nrest 
!
! parameters for lz_turb(Z) analytical below and above zm_turb turbopause
! definition of turbopause Keddy = Vumol at ~ 100-130 km 
!
   ztop = 600.e3
   zbot = 0.
   ht_bot = zm_turb/alog(lz_turb/lz_surf)
   ht_up= (ztop-zm_turb)/alog(lz_top/lz_turb)  
    
   if (jind .eq. je+1)  then
!   
! for jind = je+1 Kzu(:,je+1,:)= sum of the before calculated   kzu(is:ie,je, :) + kzw(is:ie,je, :) 
! extended corner to avoid j-haloe point update
!
    kzu(:,:) = 0.25*kzu(:,:) +0.75*kzw(:,:)
    
    return        ! extra-calculations of kzu(:,:)  for lat = je+1 are not performed
   endif  
! regular grids  js:je   
!
!  
!
! below calculations of Kzz-LES for (is:ie+1, js:je) 
! Before performing LES in Z-dir, check Convective Stability: T(z) kap2(z) 
!
   if (con_adj) call con_adjust_t2( is, ie, npz, t2, dp2, pint, kap2) 
   
   
   kzw(is:ie,:) = les_kzmin
   kzu(is:ie,:) = les_kzmin
   kzt(is:ie,:) = les_kzmin/les_prnum
   kzq(is:ie,:) = les_kzmin/les_schnum   
   kzv(is:ie+1,:) = les_kzmin  
   
    nfixed =0
    nrest=0
    
 do i=is, ie 
     ins_con(1:npz) = 0
     num_ins =0 
     kzz_edd(1:npz+1) = les_kzmin
  do k=1, npz-1
     dzg =  abs(dz2(i,k)) 
     uz =  u2(i,k)  - u2(i,k+1) 
     vz =  v2(i,k)  - v2(i,k+1) 
     tz =  (t2(i,k)  - t2(i,k+1)) /dzg 
     tc = .5*(t2(i,k)+ t2(i,k+1))
     gmid = .5*(grav2(i,k)+grav2(i,k+1))
     gcp = gmid/cp2(i,k)
     bn2 = gmid*(tz +gcp)/tc
     shr2 = max(uz*uz+vz*vz, dw2min)
     shr2 = shr2/dzg/dzg
     rik = bn2/shr2
     
  if (rik < ricrit)  then
     if(bn2 < bn2min) then
       bn_uns = bn2
       bn2 = bn2min
	if (bn_uns < 0) then 
	 ins_con(k) = 1
	 num_ins = num_ins + 1
	endif 
     endif  
  endif 
      lturbmax = dzg*0.25                 !... 1/4 of delz ....tunable       
   if (flag_lzmix == 1) then 
      zgrow = exp(zm2(i,k)/Ht2(i,k))
      zlturb=lturb*zgrow       
    else
!    
! analytical profile => lz_surf to lz_top in meters 
!  
      if( zm2(i,k) <= zm_turb) then 
        zlturb =lz_surf*exp(zm2(i,k)/ht_bot)
      else
        zlturb =lz_top*exp((ztop-zm2(i,k))/ht_up) 
      endif 
   endif
   
      zlturb = min(zlturb, lturbmax)
    
           kamp = sqrt(shr2)* zlturb*zlturb
	   ritur = bn2/shr2 *dzg *dzg
           w1 = 1./(1. + 5.*ritur)
	   w1 = w1*w1  
	   	   
!          w1 = 1./(1. + 10.*ritur)/(1.+8.*ritur)            !HB-93 of NCAR/CESM/CAM
	            
       kzz_edd(k)= min(max(kamp * w1, les_kzmin), les_kzmax)  
  enddo
       kzz_edd(npz)   =kzz_edd(npz-1) 
       kzz_edd(npz+1) =kzz_edd(npz-1)
      nfixed = nfixed + num_ins  
!
! check con_adjust_t1   nrest
!      
!    do k=2, npz
! for abs(dzg) : k- [k+1]
!       dtemp =t2(i,k)     - t2(i,k+1)
!       tz =  dtemp /dzg 
!       tc = .5*(t2(i,k)   + t2(i,k+1))  
!       gcp = grav2(i,k)/cp2(i,k)   
!       if ((tz +gcp).lt.0.) nrest= nrest+1   !convectively unstable     
!    enddo
!
! vu2, kt2, df2, kzw, kzu, kzt, kzq, kzv    
    do k=2,npz-1
        w1 = .25*( kzz_edd(k-1)+2.*kzz_edd(k)+kzz_edd(k+1) )     
        kzw(i,k) =w1
    enddo
!    
     kzw(i, npz)= .5*(kzw(i,npz-1)+kzz_edd(npz))
     kzw(i, npz+1)=kzw(i, npz)
     
    do k=2,npz 
      w1 = .5*(kzw(i,k-1)+kzw(i,k))
      kzu(i,k) =.5*(vu2(i,k-1)+ vu2(i,k)) +w1
      kzt(i,k) = w1/les_prnum + .5*(kt2(i,k-1)+ kt2(i,k))
      kzq(i,k) = w1/les_prnum + .5*(df2(i,k-1)+ df2(i,k))
    enddo   
! bottom    
    k=npz+1
     kzu(i,k) =kzu(i,k-1)
     kzt(i,k) =kzt(i,k-1)
     kzq(i,k) =kzu(i,k-1)
! top     
    k=1	 
     kzu(i,k) =kzu(i,k+1)
     kzt(i,k) =kzt(i,k+1)
     kzq(i,k) =kzu(i,k+1)  
       
 enddo         !i-loop
! 
    do i=is, ie-1
      do k=1,npz+1
       kzv(i,k) = .5*(kzw(i,k)+kzw(i+1,k)) 
      enddo  
    enddo 
    i=ie
     do k=1,npz+1
       kzv(i,k) =   kzw(i,k)
       kzv(i+1,k) = kzv(i,k)
     enddo      
!    
! kzu-updates as latitude-j-averaged
!          
   end subroutine edcoef_sles_vert   
        
!============================================================
! new configuration 2025: molecular_diffusion_run inside 
!                         molecular_diffusion.F90
!                         only horizontal diffusion
!
! LES/DAF configuration 2026: molecular and eddy 
!                         diffusion, viscosity and heat conduction
!                         inside dycore_dis_moleddy
!============================================================
  subroutine dycore_disedmol_step(ak, bk, ps, pe,  pk, pkz, cappa, u,v,w,delp,pt,q,bd,   &
                            gridstruct,flagstruct,domain,npx,npy,npz, &
                            nq,dt,it,akap,zvir,cv_air,ng,delz, grav_int, newrad, radius)
			    
  use multi_gases_mod,          only :  vicpqd, vicvqd	 
  use constants_mod,            only :  cp_air, grav 
!  use sw_core_mod,              only :	sles_tens_shear, edcoef_smag2
  use tp_core_mod, only : deln_flux_explm, deln_flux_explm_udvd_dis, deln_flux_explm_wdis  
  use tracer_manager_mod,  only: get_tracer_index  
   use field_manager_mod,  only: MODEL_ATMOS  
  implicit none 
  		    
  logical, parameter :: md_nonconsv_teold = .false.
 
! stability for subcycling:  
! unsplit vs ADIR (alternative directions)
  real,    parameter    :: stab_125 =0.125
  real,    parameter    :: stab_adir = stab_125*1.8    ! can be 0.25~0.3
  			    
  type(fv_grid_bounds_type), intent(IN) :: bd
  
  real, intent(in), dimension(npz+1):: ak, bk
    real, intent(in) :: ps  (bd%isd:bd%ied  ,bd%jsd:bd%jed)           !< Surface pressure (pascal)
    real, intent(in) :: pe  (bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  !< edge pressure (pascal +/- 2 haloe points)
    real, intent(in) :: pk  (bd%is:bd%ie,bd%js:bd%je, npz+1)          !< pe**kappa
    real, intent(in) :: pkz (bd%is:bd%ie,bd%js:bd%je,npz)             !< finite-volume mean pk   
    real, intent(in) :: cappa(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    
    
  real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(in)   :: delp(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(in)   :: delz(bd%is:,bd%js:,1:)
  
  real, intent(inout) :: w(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(inout) :: u(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz)
  real, intent(inout) :: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz)
  real, intent(inout) :: q( bd%isd:bd%ied,bd%jsd:bd%jed,npz, nq)
  
  real, intent(in) ::   grav_int(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1)
  real, intent(in) ::     newrad(bd%is:bd%ie,bd%js:bd%je,npz+1)   
  integer, intent(in) :: npx,npy,npz,nq,it,ng
  
  real, intent(in)   :: akap,zvir,dt,cv_air, radius
  
  type(fv_grid_type),  intent(INOUT), target :: gridstruct
  type(fv_flags_type), intent(IN),    target :: flagstruct
  type(domain2d), intent(inout) :: domain

  type(group_halo_update_type) :: i_pack(6)
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, npz) :: p, t, wh
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed+1, npz)  :: uh
  real, dimension (bd%isd:bd%ied+1,bd%jsd:bd%jed, npz)  :: vh
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, npz,nq) :: qh
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,npz)     :: ua, va  
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, npz) :: cond3d, diff3d, visc3d  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, npz) :: khvis, cp3d,  mol3d
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed) :: plyr, temp, visc, cond, diff, cph2, molw
!!  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed) :: daf_kh
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,nq) :: qtra

  

! Valery may 2023:
! accurate storage arrays for visous energy loss
!  
      real,  dimension(bd%isd:bd%ied,bd%jsd:bd%jed  )  :: wdis
      real,  dimension(bd%isd:bd%ied,bd%jsd:bd%jed+1  ):: udis
      real,  dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ):: vdis  
!
! LES local arrays for the horizontal mixing
!        
      real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed)     :: edvis, edthc, edmix
      
      real, dimension(npz)                             :: nsubcyle
!
! vertical dissipation => (npx, npz) 
!
      real, dimension(bd%is:bd%ie,   npz)      :: visc_heat     
      real, dimension(bd%is:bd%ie,   npz)      :: vu2, kt2, df2, kion, radnc      
      real, dimension(bd%is:bd%ie,   npz+1)    :: pint, pint_u, kzw, kzu, kzt, kzq              
      real, dimension(bd%is:bd%ie+1, npz+1)    :: pint_v, kzv  
      real, dimension(bd%is:bd%ie+1, npz+1)    :: daf_kz, p2h2e, grav2, htmod    
      real, dimension(bd%is:bd%ie,   npz)      :: pmid, pmid_u, u2, t2, w2, g2, cp2, dp2, dz2, p2h2
      real, dimension(bd%is:bd%ie,   npz)      :: u2ds
      real, dimension(bd%is:bd%ie,   npz)      :: zmid, htmid, mol2, kap2, delp2      
      real, dimension(bd%is:bd%ie+1, npz)      :: pmid_v, v2      
      real, dimension(bd%is:bd%ie, npz, nq)    :: q2     
      real, dimension(npz)                     :: akmid, bkmid    
  integer :: i, j, k, ijm
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed
  integer :: nord_md2, nord
  integer :: idir, iq, iter, ii, ind, jind
  real :: coefmax, tst, mx
  real :: dtsub
  real :: dekin
  real :: psv, psu, gmid, arad4, daf_kh, rgrav2_a, ph1, ph2, ekin_vloss 
  logical :: con_adj
  
  integer :: ind_h2o, ind_o2, ind_o3p, ind_o3, ind_he4        
  
  con_adj = .false.           ! it can be namelist-parameter ...les_conadj 
  arad4 = 4.*radius*radius
  
  rgrav2_a =Runiv/grav      ! hydrostatic dp/dr = -rhod*g_a *a2/(a+z)**2
                            ! DAF-factor a^4/(a+z)^4 for solving vertical diss-n
  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed
  ijm = (ied-isd+1)*(jed-jsd+1)
  
!------------------
!> nord = 0:   del-2     physics-based molecular and eddy dissipation/mixing
!> nord = 1:   del-4
!> nord = 2:   del-6     current for C96 higher-order numerical dissipation
!> nord = 3:   del-8 --> current for C384-768 requires more ghosting than current
!------------------   
!  
!valery : second order phys-based mixing  NORD = 0  fixed in 2024 CUA/vay
!
  nord_md2 = 0    !second order phys-based mixing  NORD = 0
  nord = nord_md2 ! no corner-updates
  
                                       call timing_on('disedmol_step')
				       
       ind_o3 = get_tracer_index (MODEL_ATMOS, 'spo3')
       ind_o2 = get_tracer_index (MODEL_ATMOS, 'spo2')
       ind_o3p = get_tracer_index (MODEL_ATMOS, 'spo')
       ind_h2o = get_tracer_index (MODEL_ATMOS, 'sphum') 
       if (md_helium) then 
         ind_he4 = get_tracer_index (MODEL_ATMOS, 'sphe') 
       else
         ind_he4 = -1 
       endif 	 
!       
!   copy arrays in 3D-local arrays that will be updated 
!   restore "pmid" from "pint" without akmid/bkmid
!      do k=1,npz
!         ph1 = ak(k  ) + bk(k  )*flagstruct%p_ref
!         ph2 = ak(k+1) + bk(k+1)*flagstruct%p_ref
!         pmid(k) = (ph2 - ph1) / log(ph2/ph1)
!      enddo

    daf_kz(:,:) = 1.0
    daf_kh =1.
    do k=1, npz
       akmid(k) = .5*(ak(k+1)+ak(k))
       bkmid(k) = .5*(bk(k+1)+bk(k))
       uh(:,:,k) = u(:,:,k)
       vh(:,:,k) = v(:,:,k)
       wh(:,:,k) = w(:,:,k)
       qh(:,:,k,:) = q(:,:,k,:)      
!
! 3D Eulerian pressure represent actual layers on which the dissipative processes work out
!       
      do j=js,je
        do i=is,ie
          t(i,j,k) = pt(i,j,k)                   ! copy kinetic temperature
!          p(i,j,k) = ps(i,j)*bkmid(k) +akmid(k)  ! mid-point pressures 
	  ph1 = ak(k)   + bk(k)  *ps(i,j)
	  ph2 = ak(k+1) + bk(k+1)*ps(i,j)
	  p(i,j,k) = (ph2 - ph1) / log(ph2/ph1)
        enddo
      enddo
    enddo

! ------- update halo to prepare for diffusion -------
                             call timing_on('COMM_TOTAL')
    call start_group_halo_update(i_pack(1), p, domain)
    call start_group_halo_update(i_pack(2), t, domain)
    call start_group_halo_update(i_pack(3), wh, domain)
    call start_group_halo_update(i_pack(4), uh, vh, domain, gridtype=DGRID_NE)
    if ( nq > 0 ) then
                                       call timing_on('COMM_TRACER')
                    call start_group_halo_update(i_pack(5), qh, domain)
                                       call timing_off('COMM_TRACER')
    endif

    call complete_group_halo_update(i_pack(1), domain)
    call complete_group_halo_update(i_pack(2), domain)
    call complete_group_halo_update(i_pack(3), domain)
    call complete_group_halo_update(i_pack(4), domain)
    if ( nq>0 ) then
                                       call timing_on('COMM_TRACER')
         call complete_group_halo_update(i_pack(5), domain)
                                       call timing_off('COMM_TRACER')
    endif
                             call timing_off('COMM_TOTAL')

    mx = 0.
!    
! Valery 2025: Comments to accelerate based on the directional split
!
! assumes that diff-stability is evaluated for 2D(X-Y) explisit diffusiion
! for "altrenative direction" algorithms  DifX => DifY:  
!                        DifmaX/Y = 0.25*[dx*dx]/dt 
! Joint/unsplit dif-step DifmaXY = 0.125*[dx*dx]/dt
!    
    coefmax= stab_adir*min(gridstruct%da_min,gridstruct%da_min_c)/abs(dt)
    
    nsubcyle(1:npz)=1
    
    do k=1, npz
!
! prepare/store molecular and eddy-LES  dissipation coefficients on the given layer
!
      do j=jsd,jed
         do i=isd,ied
            plyr(i,j) = p(i,j,k)
            temp(i,j) = t(i,j,k)
            qtra(i,j,1:nq) = qh(i,j,k,1:nq)
          enddo
      enddo
 
        
!=====================
! fill up corner for coefficient computation with values
!
! Valery: may 2025 to do check/test do we need it for NORD=0, 2-nd order diffusion 
!
      if( flagstruct%nord>0 .and. (.not. (flagstruct%regional))) then
          idir = mod(it-1,2)+1  ! alternated by 1 and 2 to avoid bias
          call copy_corners(plyr, npx, npy, idir, gridstruct%nested, bd, &
                            gridstruct%sw_corner, gridstruct%se_corner, &
                            gridstruct%nw_corner, gridstruct%ne_corner)
          call copy_corners(temp, npx, npy, idir, gridstruct%nested, bd, &
                            gridstruct%sw_corner, gridstruct%se_corner, &
                            gridstruct%nw_corner, gridstruct%ne_corner)
          do iq=1,nq
             call copy_corners(qtra(isd,jsd,iq), npx, npy, idir, &
                               gridstruct%nested, bd           , &
                            gridstruct%sw_corner, gridstruct%se_corner, &
                            gridstruct%nw_corner, gridstruct%ne_corner)
          enddo
      endif
!
! Valery: may 2025 to do check/test do we need it for NORD=0, 2-nd order diffusion 
!=====================ispo, ispo2, isphe

      call eddy_molec_coefs(ijm,ind_o3p, ind_o2, ind_he4,        &
                 plyr(isd,jsd  ),temp(isd,jsd),                  &
                 qtra(isd,jsd,1), visc(isd,jsd), cond(isd,jsd),  &
                 diff(isd,jsd), cph2(isd,jsd), molw(isd,jsd) )
		 
! option to add eddy mixing to stabilize the geo-storm dynamics, high winds > 2000 m/s
	visc3d(:,:,k)     = visc(:,:)	
	cond3d(:,:,k)     = cond(:,:)
	diff3d(:,:,k)     = diff(:,:)
	cp3d(:,:,k)       = cph2(:,:)	  
	mol3d(:,:,k)      = molw(:,:)  
      if (flagstruct%ufswam_med) then 
      
        call edcoef_sles(ng, npx, npy, bd, u(:,:,k),  v(:,:,k), w(:,:,k), delz(:,:,k), &
	                 visc, edvis, edthc, edmix, gridstruct) 
	     
        visc(:,:)     = visc(:,:) + edvis(:,:)
        cond(:,:)     = cond(:,:) + edthc(:,:)
        diff(:,:)     = diff(:,:) + edmix(:,:)
	khvis(:,:,k) = edvis(:,:)
      endif
      mx = 0.
      if(md_subcycle)then
        do j=jsd,jed
          do i=isd,ied
            tst = max(visc(i,j),cond(i,j),diff(i,j))/coefmax
            if(tst> mx) mx = tst
          enddo
        enddo
	nsubcyle(k) = mx
      endif
      
    enddo        ! loop k=1, npz
!
    if(md_subcycle) then
      call mp_reduce_max(mx)
      iter = ceiling(mx)
      if( is_master() ) write(*,*) ' Number of MD subcycles ', iter, mx, coefmax
    else
      iter = 1
    endif
!    
   do k=1, npz 
!   
! on the given layer we perform # of iter = ceiling(mx)
!     
      mx = nsubcyle(k)
    if(md_subcycle) then
      call mp_reduce_max(mx)
      iter = ceiling(mx)
      dtsub = dt/iter
!      if( is_master() ) write(*,*) ' Number of MD subcycles ', iter, mx, coefmax
    else
      iter = 1
      dtsub = dt
    endif   
           edvis(:,:) = khvis(:,:,k)
           visc(:,:) = (visc3d(:,:,k) + edvis(:,:))*dtsub
           cond(:,:) = (cond3d(:,:,k) + edvis(:,:)/les_prnum)*dtsub
           diff(:,:) = (diff3d(:,:,k) + edvis(:,:)/les_schnum)*dtsub  
!
! daf-correction in the horizontal directions
!
! to do: for D-grid u-v winds visc-coeff should be shifted:
!  kh_v(is:ie+1,:) and kh_u(:, js:je+1), including daf_kh
!	   
      if (flagstruct%ufswam_daf) then
       do j=js,je
          do i=is,ie 
	  daf_kh = arad4/(newrad(i,j,k)+newrad(i,j,k+1))**2
	  visc(i,j) =visc(i,j)*daf_kh
	  cond(i,j)= cond(i,j)*daf_kh
	  diff(i,j)= diff(i,j)*daf_kh
          enddo
        enddo	      
      endif   
           
    do ii=1,iter
      if(ii>1)then
        call complete_group_halo_update(i_pack(2), domain)
        call complete_group_halo_update(i_pack(3), domain)
        call complete_group_halo_update(i_pack(4), domain)
	
        if ( nq>0 ) then
                                           call timing_on('COMM_TRACER')
             call complete_group_halo_update(i_pack(5), domain)
                                           call timing_off('COMM_TRACER')
        endif
                                 call timing_off('COMM_TOTAL')
      endif  
!       	   
! t
        call deln_flux_explm(nord_md2,is,ie,js,je,npx,npy,cond,t(isd,jsd,k),gridstruct,bd)
!
! q should be applied only for minor species NO, O3, H2O, not for He, O, O2, N2 in the thermosphere
! 
        do iq=1,nq
!	  if ((iq == ind_h2o) .or. (iq == ind_o3) .or. (iq == ind_o3p) .or. (iq == ind_o2)) &
	  if ((iq == ind_h2o) .or. (iq == ind_o3) ) &	
          call deln_flux_explm(nord_md2,is,ie,js,je,npx,npy,diff,qh(isd,jsd,k,iq),gridstruct,bd)
        enddo

! u v 
!        call deln_flux_explm_udvd(nord_md2,is,ie,js,je,npx,npy,visc, &
!                                  uh(isd,jsd,k),vh(isd,jsd,k),gridstruct,bd)
         udis(:,:) =0.
	 vdis(:,:) =0.
        call deln_flux_explm_udvd_dis(nord_md2,is,ie,js,je,npx,npy,visc, &
                                     uh(isd,jsd,k),vh(isd,jsd,k), udis, vdis, gridstruct,bd)				  
				  
! w
!        call deln_flux_explm(nord_md2,is,ie,js,je,npx,npy,visc,wh(isd,jsd,k),gridstruct,bd)
          wdis(:,:) =0.	
        call deln_flux_explm_wdis(nord_md2,is,ie,js,je,npx,npy,visc,wh(isd,jsd,k),wdis, gridstruct,bd)	
!
! add viscous heating here.....
!  
        do j=js,je
          do i=is,ie           
            dekin = wdis(i,j) + .5*(vdis(i, j)+vdis(i+1, j))+ &
	                   .5*(udis(i,j)+udis(i,j+1))
            t(i,j,k) = t(i,j,k)+ dekin/cp3d(i,j,k)   
	            
!! cp3d=(cp_air*vicpqd(q(i,j,k,:)))
	
	  enddo
	 enddo  
	    	  
! updates       
      if(iter>1 .and. ii<iter)then
                                 call timing_on('COMM_TOTAL')
        call start_group_halo_update(i_pack(2), t, domain)
        call start_group_halo_update(i_pack(3), wh, domain)
        call start_group_halo_update(i_pack(4), uh, vh, domain, gridtype=DGRID_NE)
        if ( nq > 0 ) then
                                           call timing_on('COMM_TRACER')
                        call start_group_halo_update(i_pack(5), qh, domain)
                                           call timing_off('COMM_TRACER')
        endif
      endif
      enddo  ! iter-K dependent
!      pt(:,:,k) = t(:,:,k)
!      u(:,:,k) = uh(:,:,k)
!      v(:,:,k) = vh(:,:,k)
!      w(:,:,k) = wh(:,:,k)
!================================================================
! species to be diffused only "minor" gases
! special code for major species O-O2-N2-He, using binary diffusion
!          they form air density
!     coefficients and explicit solvers in the horizontal directions
!     d12 "O-O2" as done in GSMWAM-dycore  does not work
! 
!================================================================      
!      q(:,:,k,ind_h2o) = qh(:,:,k,ind_h2o)
!      q(:,:,k,ind_o3) =  qh(:,:,k,ind_o3)
!      q(:,:,k,ind_o3p) = qh(:,:,k,ind_o3p)  
!      q(:,:,k,ind_o2) = qh(:,:,k,ind_o2)    
!       if (md_helium) q(:,:,k, ind_he4)= qh(:,:,k,ind_he4)  
    enddo   ! k-loop Kh in horizontal directions
    
!====================================================================
!  style of the "vertical remap" fv_mapz loop and apply 
!  IMPLICIT-Zdir solvers  for molecular + LES-eddy 
!
!====================================================================

     u2ds(:,1:npz) = uh(:,1,1:npz )   
      
  do 1000 j=js,je+1
!  
! Vertical diisipation for non-staggered vars: Wh, T, qh-minor 
! wh(isd:     ,jsd:     ,1:)   !< vertical velocity (m/s)
! t(isd:ied  ,jsd:jed  ,npz)   ! kinetic temp-re
!
! pmid = p(i,j,k) pmid_u, pmid_v
! pint            pint_u, pint_v
!      
      if ( j > 1) u2ds(:,:) = u2(:,:)
      
      if ( j < je+1 ) then
       
        do i=is,ie

	 
! define interface pressure pint(is:ie, js:je, npz+1)
         do k=1, npz
	  kap2(i,k) = cappa(i,j,k)
!	  delp2(i,k) = delp(i,j,k)
	  zmid(i,k) = newrad(i,j,k) - radius
          pmid(i, k) =p(i,j,k)
	  w2(i,k) = wh(i,j,k)
	  u2(i,k) = uh(i,j,k)
	  t2(i,k) = t(i,j,k)
	  gmid = .5*(grav_int(i,j,k)+grav_int(i,j,k+1))
	  htmod(i,k) = Runiv/mol3d(i,j,k)*t2(i,k)/gmid
	  dz2(i,k) = delz(i,j,k)
	  vu2(i,k) = visc3d(i,j,k)
	  kt2(i,k) = cond3d(i,j,k)
	  df2(i,k) = diff3d(i,j,k)     
	   do iq=1, nq
	    q2(i,k,iq) = qh(i,j,k,iq)
	   enddo
	 enddo 	 
       enddo
       do i=is,ie
	 do k=1, npz+1
	  pint(i, k) = ps(i,j)*bk(k) +ak(k)
	  grav2(i,k) = grav_int(i,j,k)	  
	 enddo 
	   k=1
	   dp2(1,k) = pint(i, k+1)-pint(i, k)
	 do k=2, npz	
	  dp2(i,k) = pint(i, k+1)-pint(i, k)
	  gmid = Rgrav2_a*(t2(i,k)+t2(i,k-1))/(mol3d(i,j,k)+mol3d(i,j,k-1))
	  p2h2e(i,k) = (pint(i, k)/gmid)**2/(pmid(i,k)-pmid(i,k-1))
	 enddo 
	  k=2
	  p2h2e(i,1) = p2h2e(i,k)*(pmid(i,k)-pmid(i,k-1))/pmid(i,1)	  
        enddo	  
! define pmid_v, pint_v for is:ie+1       
        
	do i=is,ie+1
	 psv = .5*(ps(i,j)+ps(i+1,j))
	 do k=1, npz
	  pint_v(i, k) =psv*bk(k) +ak(k)
	  pmid_v(i, k) =psv*bkmid(k) +akmid(k)
	  v2(i,k) = vh(i,j,k)
	 enddo 	
	  k=npz+1
	  pint_v(i, k) =psv*bk(k) +ak(k)
	enddo 
	
      endif                              !j < je+1)
      
! define pmid_u, pint_u for is:ie and 

        do i=is,ie  
	   psu = .5*(ps(i,j)+ps(i,j+1)) 
	 do k=1, npz    
	   pint_u(i, k) =psu*bk(k) +ak(k)
	   pmid_u(i, k) =psu*bkmid(k) +akmid(k)
	   u2(i,k) = uh(i,j,k)
	 enddo 	
	  k=npz+1
	  pint_u(i, k) =psu*bk(k) +ak(k)	   	    
        enddo
	
	jind = j
!
! compute vert-diss coefficients: kzw, kzu, kzt, kzq, kzu and make "con_adj" if required
!	
        call edcoef_sles_vert(jind, je, is, ie, npz, t2, u2, v2, dp2, dz2, zmid, htmid, pmid, &
	                pint, pmid_u, pint_u, pmid_v, pint_v,cp2, grav2, vu2, kt2, df2, &
	                kap2, kzw, kzu, kzt, kzq, kzv, con_adj)	
!
     if (flagstruct%ufswam_daf) then
      do k=1, npz+1
        do i=is,ie+1 
	  daf_kz(i,k)=(radius/newrad(i,j,k))**4       ! for hydrostatic transform  'dr = -dp/rhod/gravZ'
!kzw, kzu, kzt, kzq, kzu
        enddo	
       enddo
       
      else
         daf_kz(is:ie+1,1:npz+1) =1.0     
      endif   
        call get_linear_diss(npz, is, ie, t2, zmid, kion, radnc)
	
    if (j < je+1) then 	
			
	call uvw_viscosity(dt, daf_kz,p2h2e, kion, kzw, pint,   dp2, pmid,   grav2, is, ie, isd, ied,   jsd, jed, npz, w2, 'wdis') 
	call uvw_viscosity(dt, daf_kz,p2h2e, kion, kzu, pint_u, dp2, pmid_u, grav2, is, ie, isd, ied,   jsd, jed, npz, u2, 'udis') 	 
	call uvw_viscosity(dt, daf_kz,p2h2e, kion, kzv, pint_v, dp2, pmid_v, grav2, is, ie+1, isd, ied, jsd, jed, npz, v2, 'vdis')       		
	wh(is:ie,j,1:npz)   = w2(is:ie, 1:npz)
	vh(is:ie+1,j,1:npz) = v2(is:ie+1, 1:npz)
	uh(is:ie,j,1:npz)   = u2(is:ie, 1:npz)
!	
	call tkin_conduction(dt, daf_kz, p2h2e,radnc, kzt, pint_v, dp2, pmid_v, grav2, is, ie, isd, ied, jsd, jed, npz, t2, 'tdis') 
!      		
	t(is:ie,j,1:npz) = t2(is:ie, 1:npz)
!
! molecular and eddy diffusion only for minor species, excluding the major species O-O2-N2-He
! Separate coupled binary diffusion equations are applied for "major" gases to conserve the mass 
! To do 3x3 or 2x2 coupled diffusion equations with six D(i,j) - binary coef. 
!   (o-o2, o-n2, o2-n2, o2-he, o-he, n2-he) 1D-vertical binary diffusion is in WAM physics (O-O2-N2)
!	
        do iq=1,nq
!	  if ((iq == ind_h2o) .or. (iq == ind_o3) .or. (iq == ind_o3p) .or. (iq == ind_o2)) &
	  if ((iq == ind_h2o) .or. (iq == ind_o3) ) &	
          call get_minor_diff(dt, daf_kz,p2h2e, kzq, pint, dp2, pmid,   grav2, is, ie, isd, ied, jsd, jed, npz, q2(:,:,iq), 'qdis' )
	  qh(is:ie,j,1:npz,iq) = q2(:,:,iq)
        enddo
		
   else    ! j = je+1
   
	jind = j
	u2(:,:) = uh(:,j,:)
!
! separate call only for ----> kzu(:, je+1, 1:npz)
!	
        call edcoef_sles_vert(jind, je, is, ie, npz, t2, u2, v2, dp2, dz2, zmid, htmid, pmid, &
	                pint, pmid_u, pint_u, pmid_v, pint_v,cp2, grav2, vu2, kt2, df2, &
	                kap2, kzw, kzu, kzt, kzq, kzv, con_adj)	
			
	call uvw_viscosity(dt, daf_kz, p2h2e, kion, kzu, pint_u, dp2, pmid_u, grav2, is, ie, isd, ied, jsd, jed, npz, u2, 'udis')   
	
	uh(:,j,:) = u2(:,:)
   endif 
   	
       call comp_vicous_heat(jind, je, is, ie, npz, nq, ind_o2, ind_o3p, ind_he4, &
                             u2, u2ds, v2, w2, q2, kzu, kzw, kzv, kion, dz2, visc_heat)
	 if ( j < je+1 ) then		     
           do i= is, ie
	    do k=1,npz
	     t(i,j, k) =t(i,j, k) + visc_heat(i,k)/cp3d(i,j,k)
	    enddo
	   enddo 
	  else 
           do i= is, ie
	    do k=1,npz
	     t(i,je, k) =t(i,je, k) + visc_heat(i,k)/cp3d(i,je,k)
	    enddo
	   enddo 	     	   
	  endif  
1000  continue

     do k =1, npz
      pt(:,:,k) = t(:,:,k)  
      u(:,:,k) = uh(:,:,k)
      v(:,:,k) = vh(:,:,k)
      w(:,:,k) = wh(:,:,k)
      q(:,:,k,ind_h2o) = qh(:,:,k,ind_h2o)
      q(:,:,k,ind_o3) =  qh(:,:,k,ind_o3)
!      
! major ITM-species      
!      q(:,:,k,ind_o3p) = qh(:,:,k,ind_o3p)  
!      q(:,:,k,ind_o2) =  qh(:,:,k,ind_o2)
!    if (md_helium) q(:,:,k, ind_he4)= qh(:,:,k,ind_he4)  
     enddo   
    call timing_off('disedmol_step')
    
 return
 end subroutine dycore_disedmol_step
 
 
  subroutine comp_vicous_heat(jind, je, is, ie, km, nq, ind_o2, ind_o3p, ind_he4, &
                              u2, u2ds, v2, w2, q2, kzu, kzw, kzv, kion, dz, visc_heat)
   implicit none
   integer, intent(in)   ::   jind, je, is, ie, km, nq 
   integer, intent(in)   ::   ind_o2, ind_o3p, ind_he4
   real, intent(in)      ::   u2(is:ie, km), u2ds(is:ie, km),w2(is:ie, km), v2(is:ie+1, km)   
   real, intent(in)      ::   kzu(is:ie, km+1), kzw(is:ie, km+1), kion(is:ie, km) 
   real, intent(in)      ::   kzv(is:ie, km+1),  dz(is:ie, km) 
   real, intent(in)      ::   q2(is:ie, km, nq) 
   real, intent(out)     ::   visc_heat(is:ie, km) 
   
   integer  :: i, k
   real :: wdis, vdis, udis, ek_ion, vv, uu, vvk, uuk
   real :: kwc, kuc, kvc, dzk, dwz, duz, dvz
    if (jind  < je+1) then 
     do i= is, ie
        vv = .5*(v2(i,1)+v2(i+1,1))
	uu = .5*(u2(i,1)+u2ds(i,1))
	
      do k=1, km-1 
      
        dzk = 1./abs(dz(i, k))
	
        vvk = .5*(v2(i,k+1)+v2(i+1,k+1))
	uuk = .5*(u2(i,k+1)+u2ds(i,k+1))
        ek_ion = kion(i,k)*(uu*uu + vv*vv+ w2(i,k)*w2(i,k))	
	kwc = .5*(kzw(i,k)+kzw(i,k+1))
	kvc = .25*(kzv(i,k)+kzv(i,k+1)+kzv(i+1,k)+kzv(i+1,k+1))
! later make storage kzuds for (j-1) 	
	kuc = .5*(kzu(i,k)+kzu(i,k+1))
!       kuc = .25*(kzu(i,k)+kzu(i,k+1)+kzuds(i,k)+kzuds(i,k+1) )
        dwz = (w2(i,k)-w2(i,k+1))*dzk
	wdis = kwc*dwz*dwz
	dvz = (vv-vvk)*dzk
	vdis = kvc*dvz*dvz
	duz = (uu-uuk)*dzk
	udis=  kuc*duz*duz
	
	vv = vvk
	uu = uuk
	visc_heat(i,k) = ek_ion + wdis + vdis +udis
      enddo
        k = km
      	visc_heat(i,k) =0.
     enddo 
     
      else
!      
! only U-contribution jind = je+1   visc_heat(i,k) = .5*(kion(i,k)*uu*uu + udis) 
!    
      do i= is, ie   
 	uu = .5*(u2(i,1)+u2ds(i,1))        
       do k=1, km-1 
        dzk = 1./abs(dz(i, k))       
        visc_heat(i,k) =0.
	uuk = .5*(u2(i,k+1)+u2ds(i,k+1))
	ek_ion = kion(i,k)*uu*uu
	kuc = .5*(kzu(i,k)+kzu(i,k+1))
	duz = (uu-uuk)*dzk
	udis=  kuc*duz*duz
	visc_heat(i,k) = ek_ion + udis
	uu = uuk
       enddo
        k = km
      	visc_heat(i,k) = 0.
     enddo 
     endif 
!
! 3 major gases only
!	qn2 = 1. - qo - qo2 -qhe	
!        cpx  = cpo*qo + cpo2*qo2 + cpn2*qn2 
!	if (md_helium) cpx  =	cpx + qhe * cphe
!	     
! update cp2 if it is needed after diffusion of major species ind_o2, ind_o3p, ind_he4
!   
  end subroutine comp_vicous_heat
!======================================================
! 3 subs for minor_diff, conductivity and viscosity
! 
 subroutine get_minor_diff(dt, daf, p2h2, kzz, pe, dp2, pm, grav, is, ie, isd, ied, jsd, jed, km, u, strdis)
  implicit none
  integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km  
  real, intent(in)      ::  dt 
  real, intent(inout)   ::   u(is :ie, km)
  
  real, intent(in)      ::   daf(is:ie, km+1) 
  real, intent(in)      ::   kzz(is:ie, km+1), dp2(is:ie, km), grav(is:ie, km+1) 
  real, intent(in)      ::   p2h2(is:ie, km+1)
  real, intent(in)      ::   pe(is:ie, km+1), pm(is:ie, km)
  
  character(len=*),  intent(in)    ::   strdis 
    
  integer, parameter :: dir_fac = 1
  
  integer  :: i, k
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: dta, rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km)
  
  character(len=128) :: file_dis 
  
   dta = abs(dt)
    ak(1) = 0.
    ck(km) =0.
  do i= is, ie  
    dpf = dp2(i,1)
    dpc = pe(i,1)   
    ck(1) = dta*kzz(i,2)/dpf*p2h2(i,2)*daf(i,1)
    bk(1) = ck(1) + 1. 
    
    dk(1) = u(i,1)
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)   
       ck(k) = dta/dp2(i,k)*p2h2(i,k+1)*kzz(i,k+1)  
       bk(k) = (ak(k) + ck(k))*daf(i,k) + 1.
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       bk(k) = 1. +ak(k)*daf(i,k)
    if (dir_fac == 1 ) then  
       ed(k) =  ak(k)/bk(k)
       fd(k)=   dk(k)/bk(k)
!
!   Y[k] = ed*Y[k-1] +fd[k]  ed = a(km)/b(km) fd = D(km)/b(km)
!       
       do k = km-1, 1, -1
          rinv = 1./(bk(k)-ck(k)*ed(k+1))
          ed(k) = ak(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k+1))*rinv
       enddo
!
! k=1, top-layer
!       
          u(i,1) = fd(1)
! from top to ground	  
       do k=2,km
          u(i,k) = ed(k)*u(i, k-1) + fd(k)
       enddo
    else                    !========== direction from surf to top
       ed(1) =  ck(1)/bk(1)
       fd(1)=   dk(1)/bk(1) 
       do k = 2, km, 1
          rinv = 1./(bk(k)-ak(k)*ed(k-1))
          ed(k) = ck(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k-1))*rinv
       enddo 
          u(i,km) = fd(km)     !   the bottom-layer      
       do k=km-1,1,-1
          u(i,k) = ed(k)*u(i,k+1) + fd(k)
       enddo                  
    endif    ! directional factorization: from [top=>bot] or [bot=>top]        
  enddo 
!  
 end subroutine get_minor_diff
 
 subroutine tkin_conduction(dt, daf, p2h2,kion, kzz, pe, dp2, pm, grav, is, ie, isd, ied, jsd, jed, km, u, strdis)
!
! separate work for eddy_kt x (dT/dz +grav/cp) mol_kt x dT/dz
!                    kzz=> kzz & ktmol
!
  implicit none
  integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km  
  real, intent(in)      ::  dt
  
  real, intent(inout)   ::   u(is :ie, km)
  
  real, intent(in)      ::   kion(is:ie, km),  daf(is:ie, km+1) 
  real, intent(in)      ::   kzz(is:ie, km+1), dp2(is:ie, km), grav(is:ie, km+1) 
  real, intent(in)      ::   p2h2(is:ie, km+1)
  real, intent(in)      ::   pe(is:ie, km+1), pm(is:ie, km)
  
  character(len=*),  intent(in)    ::   strdis 
    
  integer, parameter :: dir_fac = 1
  
  integer  :: i, k
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: dta, rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km)
  
  character(len=128) :: file_dis 
  
   dta = abs(dt)
    ak(1) = 0.
    ck(km) =0.
  do i= is, ie  
    dpf = dp2(i,1)
    dpc = pe(i,1)   
    ck(1) = dta*kzz(i,2)/dpf*p2h2(i,2)*daf(i,1)
    bk(1) = ck(1) + 1. +kion(i, 1)*dta
    
    dk(1) = u(i,1)
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)   
       ck(k) = dta/dp2(i,k)*p2h2(i,k+1)*kzz(i,k+1)  
       bk(k) = (ak(k) + ck(k))*daf(i,k) + 1.+kion(i, k)*dta
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       bk(k) = 1. +ak(k)*daf(i,k) +kion(i, k)*dta
    if (dir_fac == 1 ) then  
       ed(k) =  ak(k)/bk(k)
       fd(k)=   dk(k)/bk(k)
!
!   Y[k] = ed*Y[k-1] +fd[k]  ed = a(km)/b(km) fd = D(km)/b(km)
!       
       do k = km-1, 1, -1
          rinv = 1./(bk(k)-ck(k)*ed(k+1))
          ed(k) = ak(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k+1))*rinv
       enddo
!
! k=1, top-layer
!       
          u(i,1) = fd(1)
! from top to ground	  
       do k=2,km
          u(i,k) = ed(k)*u(i, k-1) + fd(k)
       enddo
    else                    !========== direction from surf to top
       ed(1) =  ck(1)/bk(1)
       fd(1)=   dk(1)/bk(1) 
       do k = 2, km, 1
          rinv = 1./(bk(k)-ak(k)*ed(k-1))
          ed(k) = ck(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k-1))*rinv
       enddo 
          u(i,km) = fd(km)     !   the bottom-layer      
       do k=km-1,1,-1
          u(i,k) = ed(k)*u(i,k+1) + fd(k)
       enddo                  
    endif    ! directional factorization: from [top=>bot] or [bot=>top]        
  enddo 
 end subroutine tkin_conduction
! 


!
 subroutine uvw_viscosity(dt, daf, p2h2, kion, kzz, pe, dp2, pm, grav, is, ie, isd, ied, jsd, jed, km, u, strdis)
 
  implicit none
  integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km  
  real, intent(in)      ::  dt
  
  real, intent(inout)   ::   u(is :ie, km)
  
  real, intent(in)      ::   kion(is:ie, km),  daf(is:ie, km+1) 
  real, intent(in)      ::   kzz(is:ie, km+1), dp2(is:ie, km), grav(is:ie, km+1) 
  real, intent(in)      ::   p2h2(is:ie, km+1)
  real, intent(in)      ::   pe(is:ie, km+1), pm(is:ie, km)
  
  character(len=*),  intent(in)    ::   strdis 
    
  integer, parameter :: dir_fac = 1
  
  integer  :: i, k
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: dta, rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km)
  logical :: mvis_debug
  character(len=128) :: file_dis 
!========================================== 
! DAF-correction for mass-based system  
!     dp_int = -rho_mid*grav_mid*delz
!
! dX    a^4
! -- = ----  x [A*X(k-1) -B*X(k) + C*X(k+1)]
! dt  (a+z)^4
!
!           wcof = .5*[rhoz(k)+rhoz(k+1)]
!           dpc= .5*[dpe(k)+dpe(k+1)]
!	    rhomol(i,k+1) =wcof*wcof*grav(i,k+1)/dpc
!===========================================  
   mvis_debug = .false.
   dta = abs(dt)
   rdt = 1./dta
   
!   if((mvis_debug) .and. (is_master())) then
!      print *, 'isd, ied, jsd, jed ', isd, ied, jsd, jed
!      print *, ' is, ie, km ',    is, ie, km 
! g/dpi [*Kz*rho*rho*g/dpm*r2/a2]
!   endif
    ak(1) = 0.
    ck(km) =0.
  do i= is, ie  
    dpf = dp2(i,1)
    dpc = pe(i,1)
!    gpf = .5*(grav(i,1)+grav(i,2))
    
    ck(1) = dta*kzz(i,2)/dpf*p2h2(i,2)*daf(i,1)
    bk(1) = ck(1) + 1. +kion(i, 1)*dta
    
    dk(1) = u(i,1)
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       gpf = kzz(i,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*p2h2(i,k+1) 
       bk(k) = (ak(k) + ck(k))*daf(i,k) + 1.+kion(i, k)*dta
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       bk(k) = 1. +ak(k)*daf(i,k) +kion(i, k)*dta
    if (dir_fac == 1 ) then  
       ed(k) =  ak(k)/bk(k)
       fd(k)=   dk(k)/bk(k)
!
!   Y[k] = ed*Y[k-1] +fd[k]  ed = a(km)/b(km) fd = D(km)/b(km)
!       
       do k = km-1, 1, -1
          rinv = 1./(bk(k)-ck(k)*ed(k+1))
          ed(k) = ak(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k+1))*rinv
       enddo
!
! k=1, top-layer
!       
          u(i,1) = fd(1)
! from top to ground	  
       do k=2,km
          u(i,k) = ed(k)*u(i, k-1) + fd(k)
       enddo
    else                    !========== direction from surf to top
       ed(1) =  ck(1)/bk(1)
       fd(1)=   dk(1)/bk(1) 
       do k = 2, km, 1
          rinv = 1./(bk(k)-ak(k)*ed(k-1))
          ed(k) = ck(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k-1))*rinv
       enddo 
          u(i,km) = fd(km)     !   the bottom-layer      
       do k=km-1,1,-1
          u(i,k) = ed(k)*u(i,k+1) + fd(k)
       enddo                  
    endif    ! directional factorization: from [top=>bot] or [bot=>top]        
  enddo
  
!    call mpp_error(FATAL,"STOP in get_molviscosity -DEBUG")                    
 end subroutine uvw_viscosity
!
!
! 
 subroutine con_adjust_t2( is, ie, km, t2, dp2, pe2, cap2)
!
! dry convective adjustment to "damp" instability in the upper layers
! optional after the vertical remap and before vert dissipation (visocity/conductivity)
! 
 
 integer            :: km, is, ie 
 real, intent(in)   ::  dp2(is:ie, km), pe2(is:ie, km+1), cap2(is:ie, km)
 real, intent(inout)::  t2(is:ie, km)
  
 integer :: i, k, kk, ktrop, jiter, niter
 
 real, dimension(km)   ::  c1dad,  c2dad, c3dad,  c4dad
 real :: zeps, rdenom, kappa, dpmid, zepsdp, zgamma, dtnext, tc1dad 
 real :: tins, atins
 logical :: stable
     ktrop = 72
     zeps = 2.0e-4
!
! make single iteration
!     
     niter= 1      
     
  do 80 i=is, ie   
       stable = .true.
       
 77    do   jiter=1, niter     
   do k= 2, ktrop 
   
         kappa =.25*(cap2(i, k)+ cap2(i,k+1))
	 dpmid = .5*(pe2(i,k+2)-pe2(i,k))
         c1dad(k) = kappa*dpmid/pe2(i,k+1)
	 c2dad(k) = (1. - c1dad(k))/(1. + c1dad(k))
	 rdenom = 1./(dp2(i,k)*c2dad(k) + dp2(i,k+1))
	 c3dad(k) = rdenom*dp2(i,k)
	 c4dad(k) = rdenom*dp2(i,k+1)
	 zepsdp = zeps*dpmid
	 
!c1dad = kappa*dpmid/pe2(k+1)
	 	 
	 zgamma = c1dad(k)*(t2(i,k) + t2(i,k+1)) 
	   tins = t2(i,k+1)-t2(i,k)
	   if (tins .ge.  (zgamma+zepsdp)) then	
	     stable = .false. 
!	     print *, ' dri_adj1 ', k, tins, t2(i,k+1), t2(i,k), zgamma, zepsdp	     
	     t2(i,k+1) = t2(i,k)*c3dad(k) + t2(i,k+1)*c4dad(k)
	     t2(i,k) =   c2dad(k)*t2(i,k+1)  
!	     atins = t2(i,k+1)-t2(i,k)
!	     print *, ' dri_adj2 ', k, atins, t2(i,k+1), t2(i,k)
	  endif  
	    
      enddo
         if(stable) goto 80  ! next i-loop
!	 zeps = zeps +zeps
!              ! next iteration to check stability
      enddo    
      
  80 continue  
  
  end subroutine con_adjust_t2   
end module dycore_dis_moleddy_mod
