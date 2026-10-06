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
      integer :: spo3, spo, spo2
      integer :: md_layers, md_tadj_layers
      real    ::  tau_visc, tau_cond, tau_diff
      real    :: md_consv_te
      logical :: md_subcycle   
!      
! 2026 md_helium flag md_helium=.true.  high-top UFSWAM configuration
!              
      logical :: md_helium
               
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
      real, parameter:: a12=9.69e18 ! O-O2 diffusion params
      real, parameter:: s12=0.774  
          
      real ::   am, fgas, a12bz, avgdbz
      real, allocatable :: visc3d(:,:,:)
!
      public molecular_diffusion_init, read_namelist_molecular_diffusion_nml
      public molecular_diffusion_coefs
      public :: visc3d, molecular_diffusion_run

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
!constants      
      a12bz = a12 * bz
      avgdbz= avgd * bz
      
#ifdef MULTI_GASES
      ind_gas_str = ind_gas
      ind_gas_end = num_gas
#else
      ind_gas_str = nwat + 1
      ind_gas_end = ncnst
#endif
      spo3  = ind_gas_str
      spo   = spo3 + 1
      spo2  = spo  + 1
     	
!	
      fgas = 1                                  ! for dry thermosphere   
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
!        d12 : diffusivity for tracers
!
!
!  Valery Yudin 2026: Optimized version for 
!                     every time step of dycore
!
!--------------------------------------------
      integer, intent(in):: dim
      
      real, intent(in):: plyr(dim), temp(dim), q(dim,*)
      
      real, intent(out):: mur(dim), lam(dim), d12(dim)
! Local:
      integer i, n
      real   am, fgas, a12bz, avgdbz
      real   qo,  qo2,  qo3,  qn2,  qh2o
      real   o_n, o2_n, o3_n, n2_n, h2o_n
      real   mu, la, t69, rho, cpx, cvx, rho1
      do n=1,dim
!check
!        if(plyr(n).le.0.0) call mpp_error(FATAL,"ERROR non positive value of plyr")
!        if(temp(n).le.0.0) call mpp_error(FATAL,"ERROR non positive value of temp")

        d12(n) = a12bz*temp(n)**s12 * temp(n)/plyr(n)
!        if( d12(n).lt.0.0 .and. is_master() ) then
!          write(*,*) 'ERROR negative d12 a12bz temp plyr s12 ', &
!                      d12(n),a12bz,temp(n),plyr(n),s12
!        endif

!        qh2o = max(0.0,q(n,1))*fgas

	
        qo   = max(0.0,q(n,spo ))
        qo2  = max(0.0,q(n,spo2))
        qo3  = max(0.0,q(n,spo3))

	
        qn2 = 1.0 - qo - qo2
        qn2  = max(min(qn2,1.0),0.0)	
!
! 3 major gases only
!		
        cpx  = cpo*qo + cpo2*qo2 + cpn2*qn2 
	 
        am = qo/amo + qo2/amo2 + qn2/amn2
        am = 1.0 / am                                   ! g/mol [29-16]
	
!	cvx  = cpx - runiv/am                           ! applied on p=const layers
	
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
  subroutine molecular_diffusion_run(u,v,w,delp,pt,pkz,cappa,q,bd,   &
                            gridstruct,flagstruct,domain,npx,npy,npz, &
                            nq,dt,it,akap,zvir,cv_air,ng,delz)
			    
  use multi_gases_mod,          only :  vicpqd, vicvqd	 
  use constants_mod,            only :  cp_air 
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
  real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(in):: pkz(bd%is:bd%ie,bd%js:bd%je,npz)
  real, intent(in):: delp(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(in):: delz(bd%is:,bd%js:,1:)
  real, intent(in):: cappa(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(inout):: w(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(inout):: u(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz)
  real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz)
  real, intent(inout):: q( bd%isd:bd%ied,bd%jsd:bd%jed,npz, nq)

  integer, intent(in):: npx,npy,npz,nq,it,ng
  
  real, intent(in)::akap,zvir,dt,cv_air
  
  type(fv_grid_type),  intent(INOUT), target :: gridstruct
  type(fv_flags_type), intent(IN),    target :: flagstruct
  type(domain2d), intent(inout) :: domain

  type(group_halo_update_type) :: i_pack(6)
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,md_layers) :: p, t, e, wh
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed+1,md_layers) :: uh
  real, dimension (bd%isd:bd%ied+1,bd%jsd:bd%jed,md_layers) :: vh
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,md_layers,nq) :: qh
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed) :: plyr, temp, visc, cond, diff
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,md_layers) :: cond3d, diff3d, visc3d_loc
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,nq) :: qtra
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed,npz) :: ua, va

! Valery may 2023:
! accurate storage arrays for visous energy loss
!  
      real,  dimension(bd%isd:bd%ied,bd%jsd:bd%jed  )  :: wdis
      real,  dimension(bd%isd:bd%ied,bd%jsd:bd%jed+1  ):: udis
      real,  dimension(bd%isd:bd%ied+1,bd%jsd:bd%jed  ):: vdis  
        
  real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed)     :: edvis, edthc, edmix
    
  integer :: i, j, k, ijm
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed
  integer :: nord_md2, nord
  integer :: idir, iq, iter, ii
  real :: coefmax, tst, mx
  real :: dtsub
  real :: dekin  
      integer :: ind_h2o, ind_o2, ind_o3p, ind_o3        
  
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
  
                                       call timing_on('d_md')
				       
       ind_o3 = get_tracer_index (MODEL_ATMOS, 'spo3')
       ind_o2 = get_tracer_index (MODEL_ATMOS, 'spo2')
       ind_o3p = get_tracer_index (MODEL_ATMOS, 'spo')
       ind_h2o = get_tracer_index (MODEL_ATMOS, 'sphum') 
       
!      t(:,:,:) = 0.
!
! copy arrays in 3D-local for the upper md-layers
!
    do k=1, md_layers
      uh(:,:,k) = u(:,:,k)
      vh(:,:,k) = v(:,:,k)
      wh(:,:,k) = w(:,:,k)
      qh(:,:,k,:) = q(:,:,k,:)
      do j=js,je
        do i=is,ie

          t(i,j,k) = pt(i,j,k)      ! kinetic temperature
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
! valery:  if md_nonconsv_teold =.true. => violation energy conservation for 2-nd order viscosity
!          conservation of energy is restored in 2023 by CUA/vay
          if( md_nonconsv_teold ) &
            e(i,j,k) = 0.5*(w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*(         &
                     u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -&
                     (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*            &
                     gridstruct%cosa_s(i,j)))

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
!=====================
      call molecular_diffusion_coefs(ijm,plyr(isd,jsd  ),temp(isd,jsd),&
                                         qtra(isd,jsd,1),visc(isd,jsd),&
                                         cond(isd,jsd  ),diff(isd,jsd))

      visc3d_loc(:,:,k) = visc(:,:)
      cond3d(:,:,k) = cond(:,:)
      diff3d(:,:,k) = diff(:,:)
!
! option to add eddy mixing to stabilize the geostorm dynamics
!      
!      call edcoef_smag2(ng, npx, npy, bd, u,  v, w, edvis, visc, gridstruct)
      
      if(md_subcycle)then
        do j=jsd,jed
          do i=isd,ied
            tst = max(visc(i,j),cond(i,j),diff(i,j))/coefmax
            if(tst> mx) mx = tst
          enddo
        enddo
      endif
      
    enddo        ! loop k=1, md_layers

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
	  if ((iq == ind_h2o) .or. (iq == ind_o3) .or. (iq == ind_o3p) .or. (iq == ind_o2)) &	
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
    enddo

    do k=1,md_layers
      pt(:,:,k) = t(:,:,k)
      u(:,:,k) = uh(:,:,k)
      v(:,:,k) = vh(:,:,k)
      w(:,:,k) = wh(:,:,k)
      q(:,:,k,ind_h2o) = qh(:,:,k,ind_h2o)
      q(:,:,k,ind_o3) = qh(:,:,k,ind_o3)
      q(:,:,k,ind_o3p) = qh(:,:,k,ind_o3p)  
      q(:,:,k,ind_o2) = qh(:,:,k,ind_o2)            
    enddo

    if( md_nonconsv_teold ) then
      do k=1,md_layers
        do j=js,je
          do i=is,ie
              e(i,j,k) = e(i,j,k) -  &
                         0.5*(w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*(         &
                         u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -&
                         (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*            &
                         gridstruct%cosa_s(i,j)))

              pt(i,j,k) = pt(i,j,k) + e(i,j,k)/cv_air    ! Valery : "bug should be cp_air with multi_gases correction
          enddo
        enddo
      enddo
    endif         !  md_nonconsv_teold =.true.

                                       call timing_off('d_md')
 return
 end subroutine molecular_diffusion_run

end module molecular_diffusion_mod
