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

!>@brief Molecular diffusion coefficients of
!>       viscosity (mur), conductivity (lam), and diffusivity (d12)
!>       and their efolding time effectiveness
!>@author H.-M. H. Juang, NOAA/NWS/NCEP/EMC
!>@Kevin Viner (NOAA/EMC) and Valery Yudin (NOAA/NESDIS/SWO-CUA) (2025) applications
!          on the Eulerian presssure layers rho*Cp*dT/dt.... with sub-cycling
!>@Valery Yudin (NOAA/NESDIS/SWO-CUA) (2026) , adding eddy LES dissipation and restoring
!          energy conservation due to the viscous loss of kinetic energy
!>@Valery Yudin (NOAA/NESDIS/SWO-CUA) (Jan 2026) adding contribution of He in the molecular
!         viscosity and heat conductivity for the top-lid (850 km) UFSWAM configurations
!>@Valery Yudin (NOAA/NESDIS/SWO-CUA) (May 2026) new dycore_dis_moleddy.F90
!                                     DAF-LES/MED configuration of eddy-mol dis-n as a sub-step
!                                     of UFSWAM-dynamical core
module molecular_diffusion_mod

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
      integer :: md_layers, md_tadj_layers
      real    :: tau_visc, tau_cond, tau_diff
      real    :: md_consv_te
      logical :: md_subcycle 
      logical :: md_helium  
      logical :: con_adj = .false.
                     
      real, parameter:: Runiv = 8314.4598 
      real, parameter:: amo=15.9994, amo2=31.9999, amo3=47.9982     !g/mol
      real, parameter::              amn2=28.013,  amh2o=18.0154    !g/mol
!< muo3 and muh2o are not precise, correct later
      real, parameter:: muo=3.9e-7, muo2=4.03e-7,  muo3=4.03e-7     !kg/m/s
      real, parameter::             mun2=3.43e-7,  muh2o=3.43e-7    !kg/m/s
!< lao3 is not precise values, but o3_n is very small
      real, parameter:: lao=75.9e-5, lao2=56.e-5,  lao3=36.e-5     !kg/m/s
      real, parameter::              lan2=56.e-5,  lah2o=55.e-5    !kg/m/s
      real, parameter:: cpo=1299.185, cpo2=918.0969, cpo3=820.2391
      real, parameter::               cpn2=1031.108, cph2o=1846.00
      real, parameter:: avgd=6.0221415e23  ! Avogadro constant
      real, parameter:: bz=1.3806505e-23   ! Boltzmann constant J/K
      real, parameter:: a12=9.69e18        ! O-O2 binary diffusion coef
      real, parameter:: s12=0.774

!
      public molecular_diffusion_init, read_namelist_molecular_diffusion_nml
      public molecular_diffusion_coefs
      public :: med_diss_step

      CONTAINS
! --------------------------------------------------------
      subroutine molecular_diffusion_init(ncnst,nwat)
!--------------------------------------------
! molecular diffusion control for each effect from namelist
! Input:
!        ncnst    : number of all prognostic tracers
!        nwat     : number of water and the end location of water
!--------------------------------------------
      integer, intent(in):: ncnst, nwat
!

      if( is_master() ) then
        write(*,*) ' molecular_diffusion is on'
        write(*,*) ' molecular_diffusion number of layers ',md_layers
        write(*,*) ' molecular_diffusion subcycling ',md_subcycle
        write(*,*) ' thermosphere adjustment number of layers ',md_tadj_layers
        write(*,*) ' energy conservation of MD (0:off, 1:on) ',md_consv_te
        write(*,*) ' viscosity    day ',tau_visc,' with effect ',tau_visc
        write(*,*) ' conductivity day ',tau_cond,' with effect ',tau_cond
        write(*,*) ' diffusivity  day ',tau_diff,' with effect ',tau_diff
      endif

#ifdef MULTI_GASES
      ind_gas_str = ind_gas
      ind_gas_end = num_gas
#else
      ind_gas_str = nwat + 1
      ind_gas_end = ncnst
#endif

      return
      end subroutine molecular_diffusion_init

! --------------------------------------------------------
      subroutine read_namelist_molecular_diffusion_nml(nml_filename,ncnst,nwat)

      character(*), intent(IN) :: nml_filename
      integer, intent(IN) :: ncnst, nwat
      integer :: ierr, f_unit, unit, ios

      namelist /molecular_diffusion_nml/ tau_visc, tau_cond, tau_diff, &
                                         md_layers, md_tadj_layers, &
                                         md_consv_te, md_subcycle, md_helium

      unit = stdlog()

      ! Read molecular_diffusion namelist
      read (input_nml_file,molecular_diffusion_nml,iostat=ios)
      ierr = check_nml_error(ios,'molecular_diffusion_nml')

      write(unit, nml=molecular_diffusion_nml)
      call molecular_diffusion_init(ncnst,nwat)

      return
      end subroutine read_namelist_molecular_diffusion_nml

! --------------------------------------------------------
      subroutine molecular_diffusion_coefs(dim,plyr,temp,q,mur,lam,d12)
!--------------------------------------------
! !OUTPUT PARAMETERS
! Input: dim : length of field
!        plyr: pressure at given layer (pascal)
!        temp: temperature (K)
!        q   : tracers, needs only gas (kg/kg)
! Ouput: mur : viscosity for momemtum
!        lam : conductivity for thermodynamics
!        d12 : effective diffusivity for tracers "O-O2"
!---------------------------------------------------------------------
! To do: extension for binary diffusion matrix of major species
!        nmsp  is a number of major species 
! dif_b  (nmsp, nmsp) is the MSP binary dif-matrix of Maxwell-Stefan Equations for "MMR"-fluxes
!        that include impact of the molecular diffusion of major species(O-O, O-O2, O-He, O-N2)
!
!        N2 = 1-O-O2-He, is not explicitly included to conserve the total mass
!
!        dmsp(nmsp, nmsp, dim) for 2 cases:
!        nmsp=2 (O-O2-N2) nmsp=3 (O-O2-He-N2), Dij -binary dif. coef.
!        molecular weight (MWA)of air and MWi of msp-gases
!-----------------------------------------------------------------------
      integer, intent(in):: dim
      real, intent(in):: plyr(dim), temp(dim), q(dim,*)
      real, intent(out):: mur(dim), lam(dim), d12(dim)
! Local:
      integer i, n, spo3, spo, spo2
      real   am, fgas, a12bz, avgdbz
      real   qo,  qo2,  qo3,  qn2,  qh2o
      real   o_n, o2_n, o3_n, n2_n, h2o_n
      real   mu, la, t69, rho, cpx, cvx, rho1

!constants
      a12bz = a12 * bz
      avgdbz= avgd * bz
      spo3  = ind_gas_str
      spo   = spo3 + 1
      spo2  = spo  + 1
      if( spo.gt.ind_gas_end ) spo=0
      if( spo2.gt.ind_gas_end ) spo2=0

      do n=1,dim
!check
        if(plyr(n).le.0.0) call mpp_error(FATAL,"ERROR non positive value of plyr")
        if(temp(n).le.0.0) call mpp_error(FATAL,"ERROR non positive value of temp")

        d12(n) = a12bz*temp(n)**s12 * temp(n)/plyr(n)
        if( d12(n).lt.0.0 .and. is_master() ) then
          write(*,*) 'ERROR negative d12 a12bz temp plyr s12 ', &
                      d12(n),a12bz,temp(n),plyr(n),s12
        endif

        fgas = 1.0
        do i=2,ind_gas_str-1
          fgas = fgas - max(0.0,q(n,i))
        enddo
        fgas = max(min(fgas,1.0),1.0E-20)         ! reasonable assured
        fgas = 1./fgas

        qh2o = max(0.0,q(n,1))*fgas
        qo   = 0.0
        qo2  = 0.0
        qo3  = 0.0
        if(spo .ne.0) qo   = max(0.0,q(n,spo ))*fgas
        if(spo2.ne.0) qo2  = max(0.0,q(n,spo2))*fgas
        if(spo3.ne.0) qo3  = max(0.0,q(n,spo3))*fgas


! reasonable values assure
!        qo   = max(min(qo ,1.0),0.0)
!        qo2  = max(min(qo2,1.0),0.0)
!        qo3  = max(min(qo3,1.0),0.0)
!        qh2o = max(min(qh2o,1.0),0.0)
	
        qn2 = 1.0 - qo - qo2
        qn2  = max(min(qn2,1.0),0.0)	
!
! 3 major gases only
!		
        cpx  = cpo*qo + cpo2*qo2 + cpn2*qn2 
	 
        am = qo/amo + qo2/amo2 + qn2/amn2
        am = 1.0 / am    ! g/mol [29-16]
	
	cvx  = cpx - runiv/am
	
        o_n   = qo   * am / amo
        o2_n  = qo2  * am / amo2
        n2_n  = qn2  * am / amn2

	
        mu = o_n*muo + o2_n*muo2 + n2_n*mun2
        la = o_n*lao + o2_n*lao2 + n2_n*lan2

        t69 = temp(n) ** 1.69                          ! T^(0.69+1)
!        rho = 1e-3 * am * plyr(n)/temp(n) / avgdbz    ! kg/m^3
        rho1 =avgdbz/plyr(n)/am
        mur(n) = mu * t69 * rho1
!vay-bug        lam(n) = la * t69 / rho / cvx
        lam(n) = la * t69 *rho1 / cpx

      enddo
     
      return
      end subroutine molecular_diffusion_coefs
      
!============================================================
! new configuration 2025: molecular_diffusion_run inside 
!  molecular_diffusion.F90
!                          only horizontal diffusion
!============================================================
  subroutine med_diss_step(u,v,w,delp,pt,pkz,cappa,q,bd,    &
                           gridstruct,flagstruct,domain,npx,npy,npz, &
                           nq,dt,it,akap,zvir,cv_air,ng, delz, ak, bk, grav_int, rad_int, pe)
			    
  use dycore_dis_moleddy_mod,   only :	eddy_molec_coefs, edcoef_sles
  use dycore_dis_moleddy_mod,   only :  comp_vicous_heat, edcoef_sles_vert
  use dycore_dis_moleddy_mod,   only :  uvw_viscosity,  tkin_conduction, get_linear_diss
  use dycore_dis_moleddy_mod,   only :  get_minor_diff, dirz_viscosity
  use multi_gases_mod,          only :  vicpqd, vicvqd	 
  use constants_mod,            only :  cp_air, grav, radius 
  
!  use sw_core_mod,              only :	sles_tens_shear, edcoef_smag2
  use tp_core_mod, only : deln_flux_explm, deln_flux_explm_udvd_dis, deln_flux_explm_wdis  
  use tracer_manager_mod,  only: get_tracer_index  
  use field_manager_mod,  only: MODEL_ATMOS  
  
  implicit none
   		    
  logical, parameter :: md_nonconsv_teold = .false.
 
! stability for subcycling:  
! using split-ADIR (alternative directions)

  real,    parameter    :: stab_125 =0.125
  real,    parameter    :: stab_adir = stab_125*1.    ! 1.8   can be 0.25~0.3
  			    
  type(fv_grid_bounds_type), intent(IN) :: bd
  type(fv_grid_type),  intent(IN), target :: gridstruct
  type(fv_flags_type), intent(IN), target :: flagstruct
  type(domain2d), intent(inout) :: domain

  type(group_halo_update_type) :: i_pack(6)
  real, intent(in)   :: pe(bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1) !< Eulerian interface pressure (pascal)
  real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(in):: pkz(bd%is:bd%ie,bd%js:bd%je,npz)
  real, intent(in):: delp(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(in):: delz(bd%is:,bd%js:,1:)
  real, intent(in):: cappa(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(inout):: w(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(inout):: u(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz)
  real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz)
  real, intent(inout):: q( bd%isd:bd%ied,bd%jsd:bd%jed,npz, nq)
  
  real, intent(in) ::   grav_int(bd%isd:bd%ied,bd%jsd:bd%jed, npz+1)
  real, intent(in) ::   rad_int (bd%isd:bd%ied,bd%jsd:bd%jed, npz+1)  
  
  integer, intent(in):: npx,npy,npz,nq,it,ng
  
  real, intent(in), dimension(npz+1):: ak, bk
  real, intent(in) :: akap,zvir,dt,cv_air

  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,md_layers) :: p, t, wh
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed+1,md_layers) :: uh
  real, dimension (bd%isd:bd%ied+1,bd%jsd:bd%jed,md_layers) :: vh
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, md_layers, nq) :: qh
!  
! 2D  
!
  real, dimension (bd%is-1:bd%ie+1,bd%js-1:bd%je+1) :: psh
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed) :: plyr, temp, visc, cond, diff
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed) :: cph2, molw 
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed) :: visc_lim, cond_lim, diff_lim 
   
! LES local arrays for the horizontal mixing
!        
  real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed)       :: edvis, edthc, edmix  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, npz) :: khvis, cp3d,  mol3d 
   
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,md_layers) :: cond3d, diff3d, visc3d_loc
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,nq) :: qtra

! Valery may 2023:
! accurate storage arrays for visous energy loss
!  
      real,  dimension(bd%isd:bd%ied,bd%jsd:bd%jed  )  :: wdis
      real,  dimension(bd%isd:bd%ied,bd%jsd:bd%jed+1  ):: udis
      real,  dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ):: vdis  
        
      real,  dimension(bd%isd:bd%ied,bd%jsd:bd%jed  )  :: wk_da


!
! Z_ADIR: vertical direction => [npx, npz] in j-loop
!
      real, dimension(bd%is:bd%ie,   npz)      :: visc_heat     
      real, dimension(bd%is:bd%ie,   npz)      :: vu2, kt2, df2, kion, radnc, kion1 
      real, dimension(bd%is:bd%ie+1, npz)      :: kionv   
      real, dimension(bd%is:bd%ie,   npz+1)    :: kzw, kzu, kzt, kzq , kzu1, kzw1
      real, dimension(bd%is:bd%ie+1, npz+1)    :: kzv , kzv1                     
      real, dimension(bd%is:bd%ie,   npz+1)    :: pint, pint_u 
                  
      real, dimension(bd%is:bd%ie+1, npz+1)    :: pint_v, grav2v,  p2h2ev, daf_kzv
         
      real, dimension(bd%is:bd%ie, npz+1)      :: daf_kz, p2h2e, grav2, p2h2eU, grav2U
      real, dimension(bd%is:bd%ie, npz, nq)    :: q2        
      real, dimension(bd%is:bd%ie,   npz)      :: pmid, pmid_u, u2, t2, w2, g2, cp2, dp2, dz2, p2h2, dp2u
      real, dimension(bd%is:bd%ie,   npz)      :: u2ds
      real, dimension(bd%is:bd%ie,   npz)      :: zmid, htmid, mol2, kap2      
      real, dimension(bd%is:bd%ie+1, npz)      :: pmid_v, v2, dp2v    
        
      real, dimension(bd%is:bd%ie,   npz+1)    :: kzu_je1, daf_kz_je1, p2h2e_je1
      real, dimension(bd%is:bd%ie,   npz)      :: kion_je1,   dp2_je1   

     
  real,  dimension(npz)                            :: akmid, bkmid, pgmid    
  real,  dimension(npz)                            :: nsubcyle
  integer, dimension(npz)                          :: isubcyle    
  
  integer :: jind
  integer :: ind_h2o, ind_o2, ind_o3p, ind_o3, ind_he4         
  integer :: i, j, k, ijm
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed
  integer :: nord_md2, nord
  integer :: idir, iq, iter, ii, md_iter1, md_iter2, md_itermax
  real :: coefmax, tst, mx, mxk
  real :: dtsub
  real :: dekin
  logical :: deep_atmd
  
  real    :: psv, psu, gmid, rgrav2_a, ph1, ph2, ekin_vloss 
  real    :: rad_mid, arad4, daf_kh
  
  arad4 = 4.*radius*radius  ! 4a^2
  rgrav2_a =Runiv/grav      ! hydrostatic dp/dr = -rhod*g_a *a2/(a+z)**2
  deep_atmd =flagstruct%deep_atmd
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
  
                                       call timing_on('med_wamrun')
				       
       ind_o3 = get_tracer_index (MODEL_ATMOS, 'spo3')
       ind_o2 = get_tracer_index (MODEL_ATMOS, 'spo2')
       ind_o3p = get_tracer_index (MODEL_ATMOS, 'spo')
       ind_h2o = get_tracer_index (MODEL_ATMOS, 'sphum') 
       
       if (md_helium) then 
         ind_he4 = get_tracer_index (MODEL_ATMOS, 'sphe') 
       else
         ind_he4 = -1 
       endif 	       
!      t(:,:,:) = 0.
!
! copy arrays in 3D-local for the upper md-layers (???) bad halo_update
!
        do j=js-1,je+1
        do i=is-1,ie+1
          psh(i,j) = pe(i, npz+1, j)
        enddo
	enddo 
      do k=1,npz
      
          akmid(k) = .5*(ak(k+1)+ak(k))
          bkmid(k) = .5*(bk(k+1)+bk(k))    
	      
         ph1 = ak(k  ) + bk(k  )*flagstruct%p_ref
         ph2 = ak(k+1) + bk(k+1)*flagstruct%p_ref
         if (ph1 .ne. 0) then 
	   pgmid(k) = (ph2 - ph1) / log(ph2/ph1) 
	 else
	   pgmid(k) = (ph2 + ph1)*.5
	 endif
!	 if( is_master() ) print *, k,  ' pgmid= ',  pgmid(k), ph1, ak(k), bk(k) 
      enddo
      
    do k=1, md_layers
      uh(:,:,k) = u(:,:,k)
      vh(:,:,k) = v(:,:,k)
      wh(:,:,k) = w(:,:,k)
      t(:,:,k)    = pt(:,:,k)      
      do iq =1,nq
        qh(:,:,k,iq) = q(:,:,k,iq)
      enddo
 
      do j=js,je
        do i=is,ie
	go to 2026
	 rad_mid = rad_int(i,j,k) + rad_int(i,j,k+1)
         wk_da(i,j)  = arad4/rad_mid/rad_mid
#ifdef MOIST_CAPPA
          p(i,j,k) = exp( log(pkz(i,j,k)) / cappa(i,j,k) )
#else
#ifdef MULTI_GASES
          p(i,j,k) = exp( log(pkz(i,j,k)) / &
                   (akap*virqd(q(i,j,k,:))/vicpqd(q(i,j,k,:))) )
#else
          p(i,j,k) = exp( log(pkz(i,j,k)) / akap )
#endif
#endif
!
! valery:  Eulerian layer ps-ak-bk
!
  2026    continue
	  ph1 = ak(k)   + bk(k)  *psh(i,j)
	  ph2 = ak(k+1) + bk(k+1)*psh(i,j)
         if (ph1 .ne. 0) then 
	   p(i,j,k) = (ph2 - ph1) / log(ph2/ph1)
	 else
	   p(i,j,k) = (ph2 + ph1)*.5
	 endif	  
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
			     
    md_iter1 = 0 
    md_iter2 = 0
    md_itermax =1 
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
    nsubcyle(1:npz)=1.
    isubcyle(1:npz)=1    
    
    do k=1,md_layers
!
! prepare molecular diffusion coefficients
!
      do j=jsd,jed
         do i=isd,ied
            plyr(i,j) = p(i,j,k)
            temp(i,j) = t(i,j,k)
            qtra(i,j,1:nq) = qh(i,j,k,1:nq)
          enddo
      enddo  
          
      do j=js,je
        do i=is,ie
	  rad_mid = rad_int(i,j,k) + rad_int(i,j,k+1)
          wk_da(i,j)  = arad4/rad_mid/rad_mid
         enddo
      enddo    	 
!=====================
! fill up corner for coefficient computation with values
!
! Valery: may 2025 we need it even for NORD=0, 2-nd order diffusion 
!
      if( .not. flagstruct%regional) then
          idir = mod(it-1,2)+1                ! alternated by 1 and 2 to avoid bias  ????
          call copy_corners(plyr, npx, npy, idir, gridstruct%nested, bd, &
                            gridstruct%sw_corner, gridstruct%se_corner, &
                            gridstruct%nw_corner, gridstruct%ne_corner)
          call copy_corners(temp, npx, npy, idir, gridstruct%nested, bd, &
                            gridstruct%sw_corner, gridstruct%se_corner, &
                            gridstruct%nw_corner, gridstruct%ne_corner)
          call copy_corners(wk_da, npx, npy, idir, gridstruct%nested, bd, &
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
!
!      call molecular_diffusion_coefs(ijm,plyr(isd,jsd  ),temp(isd,jsd),&
!                                         qtra(isd,jsd,1),visc(isd,jsd),&
!                                         cond(isd,jsd  ),diff(isd,jsd))
					 
      call eddy_molec_coefs(ijm, ind_o3p, ind_o2, ind_he4,       &
                 plyr(isd,jsd  ), temp(isd,jsd),                 &
                 qtra(isd,jsd,1), visc(isd,jsd), cond(isd,jsd),  &
                 diff(isd,jsd),   cph2(isd,jsd), molw(isd,jsd) )
		 
      visc3d_loc(:,:,k) = visc(:,:)
      cond3d(:,:,k) = cond(:,:)
      diff3d(:,:,k) = diff(:,:)
      cp3d(:,:,k)       = cph2(:,:)	  
      mol3d(:,:,k)      = molw(:,:)      
      
      if (flagstruct%ufswam_med) then 
         edvis(:,:) = 0.
         edthc(:,:) = 0.
         edmix(:,:) = 0.  
	      
        call edcoef_sles(ng, npx, npy, bd, uh(:,:,k),  vh(:,:,k), wh(:,:,k), delz(:,:,k), &
	                 visc, edvis, edthc, edmix, gridstruct) 
	     
        visc(:,:)     = visc(:,:) + edvis(:,:)
        cond(:,:)     = cond(:,:) + edthc(:,:)
        diff(:,:)     = diff(:,:) + edmix(:,:)
	khvis(:,:,k)  = edvis(:,:)
	
      endif        
!
! option to add eddy mixing to stabilize the geostorm dynamics
!      
!      call edcoef_smag2(ng, npx, npy, bd, u,  v, w, edvis, visc, gridstruct)

	  mxk  = 0.        
    if(md_subcycle)  then
        do j=jsd,jed
          do i=isd,ied	  
            tst = max(visc(i,j),cond(i,j),diff(i,j))/coefmax
            if(tst> mx)  mx = tst
	    if(tst> mxk) mxk = tst
          enddo
        enddo   
	
       isubcyle(k) =  ceiling(mxk)
       nsubcyle(k) = mxk
       if (nsubcyle(k) < 1.)      md_iter1 =md_iter1+1
       if (nsubcyle(k) .ge. 1. )  md_iter2 =md_iter2+1
      
       
      else 
!
! perform the single step without sub-cycles: visc_lim, 
!
          do j=jsd,jed
            do i=isd,ied
	    if (flagstruct%deep_atmd) then
	      daf_kh = wk_da(i,j)
	     else
	      daf_kh = 1.
	    endif
              visc_lim(i,j) = min(coefmax,visc(i,j))*dt
              cond_lim(i,j) = min(coefmax,cond(i,j))*dt
              diff_lim(i,j) = min(coefmax,diff(i,j))*dt
            enddo
          enddo 
	  
	call deln_flux_explm(nord_md2,is,ie,js,je,npx,npy, cond_lim, t(isd,jsd,k),gridstruct,bd)            
	wdis(:,:) =0.	
        call deln_flux_explm_wdis(nord_md2,is,ie,js,je,npx,npy, visc_lim, wh(isd,jsd,k),wdis, gridstruct,bd)
         udis(:,:) =0.
	 vdis(:,:) =0.
        call deln_flux_explm_udvd_dis(nord_md2,is,ie,js,je,npx,npy, visc_lim, &
                                     uh(isd,jsd,k),vh(isd,jsd,k), udis, vdis, gridstruct,bd)
				     
	iq=ind_h2o			     
	call deln_flux_explm(nord_md2,is,ie,js,je,npx,npy,  diff_lim,qtra(isd,jsd,iq),gridstruct,bd)	
        iq=ind_o3
	call deln_flux_explm(nord_md2,is,ie,js,je,npx,npy,  diff_lim,qtra(isd,jsd,iq),gridstruct,bd)
!
! add viscous heating here.....
!  
        do j=js,je
          do i=is,ie           
            dekin = wdis(i,j) + .5*(vdis(i, j)+vdis(i+1, j))+ &
	                        .5*(udis(i,j) +udis(i,j+1))
            t(i,j,k) = t(i,j,k)+ dekin/cph2(i,j)
	    qh(i,j,k,ind_h2o) = qtra(i,j,ind_h2o)
	    qh(i,j,k,ind_o3) = qtra(i,j,ind_o3)	
	  enddo
	 enddo  
	 			
      endif      ! md_subcycle
     
    enddo        ! loop k=1, md_layers
    
    
    if (.not.md_subcycle) then
                                 call timing_on('COMM_TOTAL')
        call start_group_halo_update(i_pack(2), t, domain)
        call start_group_halo_update(i_pack(3), wh, domain)
        call start_group_halo_update(i_pack(4), uh, vh, domain, gridtype=DGRID_NE)  
        call start_group_halo_update(i_pack(5), qh, domain)
	 

         call complete_group_halo_update(i_pack(2), domain)
         call complete_group_halo_update(i_pack(3), domain)
         call complete_group_halo_update(i_pack(4), domain)                                               
         call complete_group_halo_update(i_pack(5), domain) 
	                          call timing_off('COMM_TOTAL')
      goto 239 
    
    endif    
    
 
!================================================================    
! perform individual sub-sybcyles from  1-to Kmdis with 
!         limits of Kdis based on the molecular lenght path
!================================================================        
      md_itermax = maxval(isubcyle)
      if( is_master() ) write(*,*) 'md_iter1-2-max ', md_iter1, md_iter2, md_itermax
      
    if(md_subcycle) then
      call mp_reduce_max(mx)
      iter = ceiling(mx)
      if( is_master() ) write(*,*) ' Number of MD subcycles ', iter, mx, coefmax
    else
      iter = 1
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

      do k=1,  md_layers
! clean-up code take-out tau_visc tau_cond, tau_diff     
! may 2025

        if(md_subcycle)then
	   dtsub = dt/iter
          visc(:,:) = visc3d_loc(:,:,k)*dtsub
          cond(:,:) = cond3d(:,:,k)*dtsub
          diff(:,:) = diff3d(:,:,k)*dtsub
        else
!
! numerical stability restricted coeff-nts by coefmax 
! "wrong" way
!	
          do j=jsd,jed
            do i=isd,ied
              visc(i,j) = min(coefmax,visc3d_loc(i,j,k))*abs(dt)
              cond(i,j) = min(coefmax,cond3d(i,j,k))*abs(dt)
              diff(i,j) = min(coefmax,diff3d(i,j,k))*abs(dt)
            enddo
          enddo
        endif

! t
        call deln_flux_explm(nord_md2,is,ie,js,je,npx,npy,cond,t(isd,jsd,k),gridstruct,bd)

! q
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
            t(i,j,k) = t(i,j,k)+ dekin/(cp_air*vicpqd(q(i,j,k,:)))	
	  enddo
	 enddo 
	  
       	enddo    !k-loop 	      
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
      
    enddo ! iter-loop
    
239  continue             ! single time-step w/o sub-cycles for high dissipation

    do k=1,md_layers
      pt(:,:,k) = t(:,:,k)
      u(:,:,k) = uh(:,:,k)
      v(:,:,k) = vh(:,:,k)
      w(:,:,k) = wh(:,:,k)
      q(:,:,k,ind_h2o) = qh(:,:,k,ind_h2o)
      q(:,:,k,ind_o3)  = qh(:,:,k,ind_o3)
!msp      q(:,:,k,ind_o3p) = qh(:,:,k,ind_o3p)  
!msp      q(:,:,k,ind_o2) = qh(:,:,k,ind_o2)            
    enddo  

				  
  
!-------------------------------------------------------------------------------------
! Vertical portion of Z_ADIR for ufswam_med = .true. it can be performed before H_ADIR
!  to remove "dis-continuties" after the "Lagrangian_to_Eulerian" execution
!-------------------------------------------------------------------------------------

   do k=1, npz
    do i=is,ie
      u2ds(i, k) = u(i,js,k )  
    enddo
   enddo  
   
!     if( is_master() ) then  
!       write(*,*) , ' u ', maxval(u), minval(u)
!       write(*,*) , ' v ', maxval(v), minval(v)  
!       write(*,*) , ' temp ', maxval(pt), minval(pt) 
            
!     write(*,*) , ' psh ', maxval(psh), minval(psh)
!     write(*,*) , ' u2ds ', maxval(u2ds), minval(u2ds)
!     write(*,*) , ' grav_ ', maxval(grav_int(ie+1,js:je+1,:)), minval(grav_int(ie+1,js:je+1,:))
!     write(*,*) , ' rad_ ',  maxval(rad_int(ie+1,js:je+1,:)), minval(rad_int(ie+1,js:je+1,:))
!     endif
    
!$OMP parallel do default(none) shared(is,ie,js,je,npz, nq, dt, ind_o3p, ind_h2o, ind_he4, ind_o2, ind_o3, jind, con_adj, &
!$OMP                          isd, ied, jsd, jed, deep_atmd, rad_int, grav_int, pt, p, u, v, w, q, psh, ak, bk, akmid, bkmid, & 
!$OMP                         cappa, mol3d, cp3d, delz, visc3d_loc, cond3d, diff3d, u2ds, kion_je1, dp2_je1,   &
!$OMP                         rad_mid, Rgrav2_a, kzu1, kzw1, kzv1, kion1, pmid, kzu_je1, daf_kz_je1, p2h2e_je1,     &
!$OMP                         kap2, zmid, w2, u2, v2, t2, q2, gmid, htmid, dz2, vu2, kt2, df2, kion, kionv, radnc, &
!$OMP                         cp2, pint, dp2, p2h2e, psv, psu, pint_u, pint_v, pmid_u, pmid_v, dp2v, dp2u, p2h2eV, p2h2eU, &
!$OMP                         daf_kz,  daf_kzv, grav2, grav2u, grav2v, kzw, kzu, kzt, kzq, kzv, visc_heat)
      
  do 1000 j=js,je+1
!  
! 
! pmid = p(i,j,k) pmid_u, pmid_v
! pint            pint_u, pint_v
!      
      if ( j > js) u2ds(is:ie,1:npz) = u2(is:ie,1:npz)
      
      if ( j < (je+1) ) then
       
        do i=is,ie
         do k=1, npz
	  kap2(i,k) = cappa(i,j,k)
	  zmid(i,k) = rad_int(i,j,k) - radius
          pmid(i, k) = p(i,j,k)
	  w2(i,k) = w(i,j,k)
	  t2(i,k) = pt(i,j,k)
	  gmid = .5*(grav_int(i,j,k)+grav_int(i,j,k+1))
	  htmid(i,k) = Runiv/mol3d(i,j,k)*t2(i,k)/gmid
	  dz2(i,k) = -delz(i,j,k)
	  vu2(i,k) = visc3d_loc(i,j,k)
	  kt2(i,k) = cond3d(i,j,k)
	  df2(i,k) = diff3d(i,j,k) 
	  cp2(i,k) =cp3d(i,j,k)    
	   do iq=1, nq
	    q2(i,k,iq) = q(i,j,k,iq)
	   enddo
	 enddo 	 
       enddo
! cell-centered   is:ie    
       do i=is,ie
	 do k=1, npz+1
	  pint(i, k) = psh(i,j)*bk(k) +ak(k)
	  grav2(i,k) = grav_int(i,j,k)	  
	 enddo 
	   k=1
	   dp2(i,k) = pint(i, k+1)-pint(i, k)
	 do k=2, npz	
	   dp2(i,k) = pint(i, k+1)-pint(i, k)
	  gmid = Rgrav2_a*(t2(i,k)+t2(i,k-1))/(mol3d(i,j,k)+mol3d(i,j,k-1))   !RT/(g*mu)
	  p2h2e(i,k) = (pint(i, k)/gmid)**2/(pmid(i,k)-pmid(i,k-1))
	 enddo 
	  k=2
	  p2h2e(i,1) = p2h2e(i,k)*(pmid(i,k)-pmid(i,k-1))/pmid(i,2)
	  p2h2e(i,npz+1) = p2h2e(i,npz)  
        enddo
		  
! define pmid_v, pint_v for is:ie+1     
  
         call get_linear_diss(npz, is, ie, t2, zmid, kion, kionv, radnc)
	 
	do i=is,ie+1
	  psv = .5*(psh(i-1,j)+psh(i,j))
	  pint_v(i, npz+1) =psv	 
	 do k=1, npz
	  pint_v(i, k) =psv*bk(k) +ak(k)
	  pmid_v(i, k) =psv*bkmid(k) +akmid(k)
	  v2(i,k) = v(i,j,k)
	  grav2v(i,k) = .5*(grav_int(i-1,j,k)+grav_int(i,j,k))  
	 enddo
	 	   k=1
	   dp2v(i,k) = pint_v(i, k+1)-pint_v(i, k)
	 do k=2, npz	 
	  if (i.ne.ie+1) then 
	    gmid = Rgrav2_a*(t2(i,k)+t2(i,k-1))/(mol3d(i,j,k)+mol3d(i,j,k-1))
	    else 
	    gmid = Rgrav2_a*(t2(i-1,k)+t2(i-1,k-1))/(mol3d(i,j,k)+mol3d(i,j,k-1))
	  endif  
	  p2h2ev(i,k) = (pint_v(i, k)/gmid)**2/(pmid_v(i,k)-pmid_v(i,k-1))
	  dp2v(i,k) = pint_v(i, k+1)-pint_v(i, k)	
	 enddo 	
	  k=npz+1
	  grav2v(i,k) = .5*(grav_int(i,j,k)+grav_int(i-1,j,k))	 
	  k=2
	  p2h2eV(i,1) = p2h2eV(i,k)*(pmid_v(i,k)-pmid_v(i,k-1))/pmid_v(i,2)
	  p2h2eV(i,npz+1) = p2h2eV(i,npz)  	
	enddo 
	
      endif                              !j < je+1)
      
! define pmid_u, pint_u for is:ie and    js:je+1

        do i=is,ie  
	   psu = .5*(psh(i,j)+psh(i,j-1)) 
	  do k=1, npz    
	   pint_u(i, k) =psu*bk(k) +ak(k)
	   pmid_u(i, k) =psu*bkmid(k) +akmid(k)
	   u2(i,k) = u(i,j,k)
	   grav2u(i,k) = .5*(grav_int(i,j-1,k)+grav_int(i,j,k))
	  enddo 
	  pint_u(i, npz+1) =psu	  
	 do k=2, npz	 
	  gmid = Rgrav2_a*(t2(i,k)+t2(i,k-1))/(mol3d(i,j,k)+mol3d(i,j,k-1))
	  p2h2eU(i,k) = (pint_u(i, k)/gmid)**2/(pmid_u(i,k)-pmid_u(i,k-1))	
	 enddo 	
	  k=npz+1
	  grav2u(i,k) = .5*(grav_int(i,j,k)+grav_int(i,j-1,k))	 
	  k=2
	  p2h2eU(i,1) = p2h2eU(i,k)*(pmid_U(i,k)-pmid_U(i,k-1))/pmid_U(i,2)
	  p2h2eU(i,npz+1) = p2h2eU(i,npz)  	 	   	    
        enddo
	
	jind = j

     if (deep_atmd) then
      do k=1, npz
        do i=is,ie 
	  rad_mid = .5*(rad_int(i,j,k) + rad_int(i,j,k+1))
	  daf_kz(i,k)=(radius/rad_mid)**4                   ! for hydrostatic transform  'dr = -dp/rhod/gravZ'
	  daf_kzv(i,k)=(radius/rad_mid)**4 	  
!kzw, kzu, kzt, kzq, kzu
        enddo
	   i= ie+1
!	   rad_mid =rad_int(i,j,k)       !(radius/rad_mid)**4 rad_midV can be corrected by geometry...larger
	   daf_kzv(i,k)=daf_kzv(i-1,k)   ! 
       enddo   
! surface         	
	daf_kz(ie:is,npz+1)=1.
	daf_kzv(ie:is+1,npz+1)=1.      
      else
         daf_kz(is:ie,1:npz+1)    =1.0 
	 daf_kzv(is:ie+1,1:npz+1) =1.0     
      endif   
	
    if (j < (je+1) ) then 	
!
! compute vert-diss coefficients: kzw, kzu, kzt, kzq, kzu and make "con_adj" if required
!	
        call edcoef_sles_vert(jind, je, is, ie, npz, t2, u2, v2, dp2, dz2, zmid, htmid, pmid, &
	                pint, pmid_u, pint_u, pmid_v, pint_v,cp2, grav2, vu2, kt2, df2, &
	                kap2, kzw, kzu, kzt, kzq, kzv, con_adj)	
!			
! above the tropopause pint/pmid/dp2(x,y) = const, changes only due to T(x,y), Molw(x,y), grav(x,y)
!

!	if( is_master() .and. j == je ) then 
!	write(*,*) , ' u2 ', maxval(u2), minval(u2)
!	write(*,*) , ' w2 ', maxval(w2), minval(w2)	
!	write(*,*) , ' v2 ', maxval(v2), minval(v2)
!	write(*,*) , ' t2 ', maxval(t2), minval(t2)
!	write(*,*) , ' kion ', maxval(kion), minval(kionv)
!	write(*,*) , ' vu2 ', maxval(vu2), minval(vu2), maxval(kt2), minval(df2)
!	write(*,*) , ' kzw ', maxval(kzw),  maxval(kzt),  maxval(kzu), maxval(kzv)
!	write(*,*) , ' daf_kz ', maxval(daf_kz), minval(daf_kz), maxval(daf_kzv), minval(daf_kzv)
	
!	write(*,*) , ' pmid ',maxval(pmid), minval(pmid), maxval(pmid_v), minval(pmid_v)
!	write(*,*) , ' pint ',maxval(pint), minval(pint), maxval(pint_v), minval(pint_v)
!	write(*,*) , ' dp2  ',maxval(dp2),  minval(dp2)
!	write(*,*) , ' p2h2e  ',  maxval(p2h2e), minval(p2h2e)
!	endif	

				
        call dirz_viscosity(dt, daf_kz, p2h2e, kion, kzw, dp2, is, ie, npz, w2, 'wdis')
	call dirz_viscosity(dt, daf_kz, p2h2e, kion, kzu, dp2, is, ie, npz, u2, 'udis')
	call dirz_viscosity(dt, daf_kzv,p2h2ev, kionv, kzv, dp2v, is, ie+1, npz, v2, 'vdis')
		
!	if( is_master() .and. j== je ) then 
!	write(*,*) , 'af u2 ', maxval(u2), minval(u2)
!	write(*,*) , 'af w2 ', maxval(w2), minval(w2)		
!	write(*,*) , 'af v2 ', maxval(v2), minval(v2)
!	endif
	
	w(is:ie,j,1:npz)   = w2(is:ie, 1:npz)
	v(is:ie+1,j,1:npz) = v2(is:ie+1, 1:npz)
	u(is:ie,j,1:npz)   = u2(is:ie, 1:npz)
!	
!	call tkin_conduction(dt, daf_kz, p2h2e,radnc, kzt, pint_v, dp2, pmid_v, grav2, is, ie, isd, ied, jsd, jed, npz, t2, 'tdis')
 
	call dirz_viscosity(dt, daf_kz, p2h2e, radnc, kzt, dp2, is, ie, npz, t2, 'tdis')
		
!	if( is_master() .and. j==je ) then write(*,*) , 'af t2 ', maxval(t2), minval(t2)			

!      		
	pt(is:ie,j,1:npz) = t2(is:ie, 1:npz)
	
! TO DO
! molecular and eddy diffusion only for minor species, excluding the major species O-O2-N2-He
! Separate coupled binary diffusion equations are applied for "major" gases to conserve the mass 
! To do 3x3 or 2x2 coupled diffusion equations with six D(i,j) - binary coef. 
!   (o-o2, o-n2, o2-n2, o2-he, o-he, n2-he) 1D-vertical binary diffusion is in WAM physics (O-O2-N2)
!	
        do iq=1,nq
	  if ((iq == ind_h2o) .or. (iq == ind_o3) ) &	
          call get_minor_diff(dt, daf_kz, p2h2e, kzq, pint, dp2, pmid,   grav2, is, ie, &
	                      isd, ied, jsd, jed, npz, q2(:,:,iq), 'qdis' )
	  do k=1, npz		      
	  q(is:ie,j,k,iq) = q2(is:ie,k,iq)
	  enddo
        enddo
	if (jind == js) then
	   kzu1 = kzu
	   kzv1 = kzv
	   kzw1 =kzw
	   kion1=kion	   
	endif
	if (jind == je) then	
	  kzu_je1 = kzu
	  dp2_je1 = dp2
	  p2h2e_je1 = p2h2e
	  daf_kz_je1 = daf_kz 
	  kion_je1 = kion
	endif	
		
   else    ! j = je+1 last column for u2 on D-grid
   
	jind = j
	u2(is:ie, 1:npz) = u(is:ie,j,1:npz)
!
! separate call only for ----> kzu(:, je+1, 1:npz)
!	
!        call edcoef_sles_vert(jind, je, is, ie, npz, t2, u2, v2, dp2, dz2, zmid, htmid, pmid, &
!	                pint, pmid_u, pint_u, pmid_v, pint_v,cp2, grav2, vu2, kt2, df2, &
!	                kap2, kzw, kzu, kzt, kzq, kzv, con_adj)		
	
!	call uvw_viscosity(dt, daf_kz, p2h2e, kion, kzu, pint_u, dp2, pmid_u, grav2, is, ie, isd, ied, jsd, jed, npz, u2, 'udis')   
        call dirz_viscosity(dt, daf_kz_je1, p2h2e_je1, kion_je1, kzu_je1, dp2_je1, is, ie, npz, u2, 'udis')
           do i= is, ie
	    do k=1,npz		
	     u(i,j,k) = u2(i,k)
	    enddo
	  enddo
   endif 
!   points: j=1 and j= je+1 are special; for  j=2:je uc=.5*(u2ds(j-1)+ u2(j))
!           for j=1 no impact of u2(j=2)
!               j=je+1
!  comp_vicous_heat: u(je+1,:,:) contribution should be added, it needs separate js:je loop
!
   	
      if ( j < je .and. j > js ) then
       call comp_vicous_heat(jind, je, is, ie, npz, nq, ind_o2, ind_o3p, ind_he4, &
                             u2, u2ds, v2, w2, q2, kzu, kzw, kzv, kion, dz2, visc_heat)	 		     
           do i= is, ie
	    do k=1,npz
	     pt(i,j, k) =pt(i,j, k) + visc_heat(i,k)/cp2(i,k)
	    enddo
	   enddo 
!              add contribution of u(:, je+1,:) in viscous heat	   	   
      else 
	   if (j == je+1 ) then 
            call comp_vicous_heat(jind, je, is, ie, npz, nq, ind_o2, ind_o3p, ind_he4, &
                             u2, u2ds, v2, w2, q2, kzu, kzw, kzv, kion, dz2, visc_heat)	   
           do i= is, ie
	    do k=1,npz
!	     pt(i,je, k) =pt(i,je, k) + visc_heat(i,k)/cp3d(i,je, k)
	    enddo
	   enddo 
	endif	     	   
     endif  
     
1000  continue 

            j = js
	do  k=1, npz
	 do i=is,ie 
	  u2(i,k) = u(i,j+1,k)
	  u2ds(i,k) = u(i,j,k)
	  v2(i,k) = v(i,j,k)
	  w2(i,k) = w(i,j,k)
	  do iq =1, nq
	   q2(i,k,iq) = q(i,j,k,iq)
	  enddo
	  dz2(i,k) = -delz(i,j,k)
	enddo 
	enddo 	  
       call comp_vicous_heat(js, je, is, ie, npz, nq, ind_o2, ind_o3p, ind_he4, &
                             u2, u2ds, v2, w2, q2, kzu1, kzw1, kzv1, kion1, dz2, visc_heat)
			     
           do i= is, ie
	    do k=1,npz
	     pt(i,j, k) =pt(i,j, k) + visc_heat(i,k)/cp3d(i,j,k)
	    enddo
	   enddo  

            call timing_off('med_wamrun')
 return
 end subroutine med_diss_step

end module molecular_diffusion_mod
