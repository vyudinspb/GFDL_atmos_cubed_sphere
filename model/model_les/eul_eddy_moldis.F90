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

!>@brief Molecular and eddy dissipation
!>       viscosity, heat conductivity, and diffusion 
!>       horizontal and vertical directions
!>       call from fv_dynmaics  after remap ans kap_cor
!>
!>  call eul_eddy_moldis(u, v, w, delp, pt, pkz, cappa, q, bd,   &
!>            gridstruct, flagstruct, domain, i_pack, npx, npy, npz, nq, mdt, n_map, last_step, akap, zvir, cv_air)
!>   
!> Valery Yudin (2023) and Henry Juang (2020, horiz mo-diff)

module eul_eddy_moldis_mod

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


use constants_mod,       only: kappa, pi=>pi_8, omega, rdgas, grav, rvgas, cp_air, radius

  public eul_eddy_moldis, get_coef_mdif
  
  subroutine eul_eddy_moldis(u, v, w, delp, delz, pt, q,   pkz,cappa, pe, ps, phis,    &
                            bd, gridstruct,flagstruct,domain,   i_pack,npx,npy,npz, &
                            nq, mdt,it, last_step, akap,zvir,cv_air)
			    
  integer, intent(in):: npx,npy,npz,nq,it
  			    
  type(fv_grid_bounds_type), intent(IN) :: bd
  real, intent(in)   ::  mdt
  real, intent(in)   ::  pe(bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1) 
  
  real, intent(in)   ::  ps(bd%isd:bd%ied  ,bd%jsd:bd%jed)  
  real, intent(in)   ::  phis(bd%isd:bd%ied,bd%jsd:bd%jed)       !< Surface geopotential (g*Z_surf)  
  
  real, intent(inout):: pkz(bd%is:bd%ie,bd%js:bd%je,npz)
  real, intent(inout):: pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  
  real, intent(in)   :: delp(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(in)   :: delz(bd%isd:bd%ied,bd%jsd:bd%jed,npz)  
  
  real, intent(inout):: cappa(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, intent(inout):: w(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  
  real, intent(inout):: u(bd%isd:bd%ied,bd%jsd:bd%jed+1,npz)
  real, intent(inout):: v(bd%isd:bd%ied+1,bd%jsd:bd%jed,npz)
  
  real, intent(inout):: q( bd%isd:bd%ied,bd%jsd:bd%jed,npz, nq)
  
  logical, intent(in):: last_step

  real, intent(in):: akap, zvir, cv_air
  
  type(fv_grid_type),  intent(INOUT), target :: gridstruct
  type(fv_flags_type), intent(IN),    target :: flagstruct
  
  type(domain2d), intent(inout) :: domain
  type(group_halo_update_type), intent(inout) :: i_pack(*)

  real :: pkzf(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed) :: p, t, dekin
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed) :: wdis
  real ::  vdis(bd%isd:bd%ied+1,bd%jsd:bd%jed), udis(bd%isd:bd%ied,bd%jsd:bd%jed+1) 
  real :: wrk2, wrk1    
  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, km+1) ::  vu3d
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, km+1) ::  kt3d 
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, km+1) ::  df3d  
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, km+1) ::  ktur3d   
  real, dimension (bd%isd:bd%ied,bd%jsd:bd%jed, km+1) ::  cp_mu3d
       
  real, dimension(is:ie+1, km) :: dpv, dpu   
  real, dimension(is:ie,km+1)  :: vumol, ktmol, kteddy, dfmol, rhomol, wgrav, cp_mu
  real, dimension(is:ie,km+1)  :: zgeo   
  real, dimension(isd:ied, jsd:ed) :: vu2, kt2, ktur2, df2, cp2   
  integer :: i, j, k
  integer :: is,  ie,  js,  je
  integer :: isd, ied, jsd, jed
  integer :: ind_h2o, ind_o2, ind_o3p, ind_o3   
  real :: pe1(bd%is : bd%ie, km+1)
  real :: dpe(bd%is : bd%ie, km)
  
  real :: gravr2, pdp
    
      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed
      
       ind_o3 = get_tracer_index (MODEL_ATMOS, 'spo3')
       ind_o2 = get_tracer_index (MODEL_ATMOS, 'spo2')
       ind_o3p = get_tracer_index (MODEL_ATMOS, 'spo')
       ind_h2o = get_tracer_index (MODEL_ATMOS, 'sphum') 
       
                                       call timing_on('d_md')
! -----------------------------------------------------
! ------- update halo to prepare for diffusion -------
    pkzf = 0.0
    do k=1,npz
       do j=js,je
        do i=is,ie
          if (pkz(i,j,k) <= 0.) print *, ' mdif_run pkz < 0  see fv_mapz' 
          pkzf(i,j,k) = pkz(i,j,k)
	 enddo
       enddo
    enddo

                             call timing_on('COMM_TOTAL')
    call start_group_halo_update(i_pack(1),delp,  domain, complete=.false.)
    call start_group_halo_update(i_pack(1), pt,   domain, complete=.true.)
    call start_group_halo_update(i_pack(2),pkzf,  domain)
    call start_group_halo_update(i_pack(7), w, domain)
    call start_group_halo_update(i_pack(8), u, v, domain, gridtype=DGRID_NE)
    if ( nq > 0 ) then
                                       call timing_on('COMM_TRACER')
                    call start_group_halo_update(i_pack(10), q, domain)
                                       call timing_off('COMM_TRACER')
    endif
#ifdef MOIST_CAPPA
    call start_group_halo_update(i_pack(12), cappa, domain)
#endif

    call complete_group_halo_update(i_pack(1), domain)  ! delp. pt
    call complete_group_halo_update(i_pack(2), domain)  ! pkzf
    call complete_group_halo_update(i_pack(7), domain)  ! w
    call complete_group_halo_update(i_pack(8), domain)
    if ( nq>0 ) then
                                       call timing_on('COMM_TRACER')
         call complete_group_halo_update(i_pack(10), domain)
                                       call timing_off('COMM_TRACER')
    endif
#ifdef MOIST_CAPPA
    call complete_group_halo_update(i_pack(12), domain)
#endif
                             call timing_off('COMM_TOTAL')

    if( flagstruct%nord>0 .and. (.not. (flagstruct%regional))) then
        i=mod(it-1,2)+1        ! alternatively to avoid bias
        do k=1,npz
        call copy_corners(pt(isd,jsd,k), npx, npy, i, gridstruct%nested, bd, &
                          gridstruct%sw_corner, gridstruct%se_corner, &
                          gridstruct%nw_corner, gridstruct%ne_corner)
        call copy_corners(pkzf(isd,jsd,k), npx, npy, i, gridstruct%nested, bd, &
                          gridstruct%sw_corner, gridstruct%se_corner, &
                          gridstruct%nw_corner, gridstruct%ne_corner)
        call copy_corners(cappa(isd,jsd,k), npx, npy, i, gridstruct%nested,bd, &
                          gridstruct%sw_corner, gridstruct%se_corner, &
                          gridstruct%nw_corner, gridstruct%ne_corner)
        enddo
    endif

      gravr2 = grav *radius *radius
      
     do 1000 j=js,je+1
        do i=is, ie
	
	 zgeo(i, km+1) = phis(i,j)
	 wgrav(i,km+1) = grav      ! 9.8065
	 
	 do k=km, 1, -1
	  zgeo(i, k) = zgeo(i, k+1) -delz(i,j,k)
	  
	  zgeo(i, k) = zgeo(i, k+1) + rdmulti * pdp * pt(i,j,k)
	  
	  wgrav(i,k) = gravR2/(zgeo(i, k)+ radius)**2
	enddo
	
      do k=1,km+1
        do i=is,ie
	   
           pe1(i,k) = pe(i,k,j)
	   dpe(i,k) =delp(i,j,k)
        enddo
      enddo
       if ( mdz_tk .and. (j < je+1) ) then
         call  coef_edmoldiff(mdt, vumol, ktmol, kteddy, dfmol, rhomol, wgrav,   &
	            q, pt, w, u, v, dp2, pe2, dpu(is:ie,:), dpv, pe3u(is:ie,:), pe3, cp_mu,   &
                              j, je, is, ie, isd, ied, jsd, jed, km, nq,  &
	                      ind_h2o, ind_o2, ind_o3p, ind_o3 ) 
			      
	 vu3d(is:ie, j, 1:km+1) = vumol
	 kt3d(is:ie, j, 1:km+1) = ktmol
	 df3d(is:ie, j, 1:km+1) = dfmol
	 ktur3d(is:ie, j, 1:km+1) = kteddy
	 cp_mu3d(is:ie, j, 1:km+1) = cp_mu	 		      
         else 
	 
        call  coef_edmoldiff(mdt, vumol, ktmol, kteddy, dfmol, rhomol, wgrav,   &
	  q, pt, w, u, v, dpU(is:ie,:), pe3U(is:ie,:), dpu(is:ie,:), dpv, pe3u(is:ie,:), pe3, cp_mu,   &
                              j, je, is, ie, isd, ied, jsd, jed, km, nq,  &
	                      ind_h2o, ind_o2, ind_o3p, ind_o3 )
			      
	 vu3d(is:ie, j, 1:km+1) = vumol
	 		 			      	 			      
       endif      
     enddo


!$OMP parallel do default(none) shared(npz,flagstruct,gridstruct,bd,      &
!$OMP                       last_step, it,mdt,is,ie,js,je,isd,ied,jsd,jed, &
!$OMP                                  pt,u,v,w,q,pkz,pkzf,cappa,akap,nq, &
!$OMP                            wrk2, wrk1 ,zvir,cv_air,md_layers,md_consv_te) &
!$OMP                          private(k,i,j,p,t,e, dekin, wdis, vdis, udis)
! ----------------
    do k=1, md_layers
! ----------------

! ------- prepare p and t for molecular diffusion coefficients

       do j=jsd,jed
          do i=isd,ied
           if (last_step) then 
	       t(i,j) = pt(i,j,k)
	   else
	       t(i,j) = pt(i,j,k) * pkzf(i,j,k)
	   endif 
	     
	  if (.not.last_step) then  
#ifdef MULTI_GASES
             t(i,j) = t(i,j) / virq(q(i,j,k,:))
#else
             t(i,j) = t(i,j) / (1+zvir*q(i,j,k,1))
#endif
            endif
!   zvir =  rvgas/rdgas - 1.     !< = 0.607789855
	    
#ifdef MOIST_CAPPA
#ifdef MULTI_GASES

             p(i,j) = exp( log(pkzf(i,j,k)) / &
                      (akap*virqd(q(i,j,k,:))/vicpqd(q(i,j,k,:))) )
#else
             p(i,j) = exp( log(pkzf(i,j,k)) / cappa(i,j,k) )
#endif	     		      
#else
             p(i,j) = exp( log(pkzf(i,j,k)) / akap )
#endif

	    if (t(i,j) .le. 50.) then
	      wrk1=virq(q(i,j,k,:))
	      wrk2=   akap*virqd(q(i,j,k,:))/vicpqd(q(i,j,k,:))	      
	      print *, ' mdif_run t<50K k=', k, p(i,j), pt(i,j,k), pkzf(i,j,k)
	      print *, ' mdif_run t<50K wrks', wrk1 , wrk2
	    endif
	    
          enddo
        enddo

! compute molecular diffusion with implicit time and dimensional splits

           vu2(is:ie, js:je) = vu3d(is:ie, js:je, k)
	   kt2(is:ie, js:je) = kt3d(is:ie, js:je, k)
	   ktur2(is:ie, js:je) = ktur3d(is:ie, js:je, k)
	   df2(is:ie, js:je) = df3d(is:ie, js:je, k)
	   cp2(is:ie, js:je)  = cp3d(is:ie, js:je, k)
	   
       call edmoldiff_xy( p(isd,jsd), t(isd,jsd),         &
                  u(isd,jsd,k), v(isd,jsd,k), w(isd:,jsd:,k), q,               &
		  vu2, kt2, ktur2, df2, cp2, wdis(isd,jsd), vdis(isd,jsd), udis(isd,jsd),&
                  it, nq, k, npz, mdt, gridstruct, flagstruct, bd)

       do j=js,je
          do i=is,ie
	  	     
      if (.not.last_step) then  
#ifdef MULTI_GASES
             t(i,j) = t(i,j) * virq(q(i,j,k,:))
#else
             t(i,j) = t(i,j) * (1+zvir*q(i,j,k,1))
#endif
            pt(i,j,k) = t(i,j) / pkzf(i,j,k)
	 else  
            pt(i,j,k) = t(i,j)	
	 endif    
        enddo
       enddo

! -------------------------------------------------
    enddo       ! k loop of 2d molecular diffusion
! -------------------------------------------------
                                       call timing_off('d_md')
 return
 end subroutine   eul_eddy_moldis
!=============================== 
!
 subroutine get_coef_mdif(is, ie, km, t2, qo, qo2, qo3, qh2o, dp2, pe, grav,  &
                          vumol, ktmol, dfmol, rhomol,cp_mu)	
 use constants_mod,       only:  rdgas
 implicit none
 
 integer, intent(in)   ::  is, ie , km
 real, dimension(is:ie,km) :: qo, qo2, qh2o, qo3  
 real, dimension(is:ie,km) :: dp2, t2
 real, dimension(is:ie,km+1) :: pe, grav,  vumol, ktmol, dfmol, rhomol
 real, dimension(is:ie,km+1) :: cp_mu
 real, parameter::  amo=15.9994, amo2=2.*amo, amo3= 3.*amo      
 real, parameter::  amn2=28.013,  amh2o=18.0154    !g/mol
 
!< muo3 and muh2o are not precise, correct later
 real, parameter:: muo=3.9e-7, muo2=4.03e-7,  muo3=4.03e-7     !kg/m/s
 real, parameter::             mun2=3.43e-7,  muh2o=3.43e-7    !kg/m/s
!< lao3 is not precise values, but o3_n is very small
  real, parameter:: lao=75.9e-5, lao2=56.e-5,  lao3=36.e-5     !kg/m/s
  real, parameter::              lan2=56.e-5,  lah2o=55.e-5    !kg/m/s
  real, parameter:: cpo=1299.185, cpo2=918.0969, cpo3=820.2391
  real, parameter:: cpn2=1031.108, cph2o=1846.00
  real, parameter:: avgd=6.0221415e23  ! Avogadro constant
  real, parameter:: bz=1.3806505e-23   ! Boltzmann constant J/K
  real, parameter:: a12=9.69e18 ! O-O2 diffusion params
  real, parameter:: s12=0.774, s121 = s12+1.
  real, parameter:: ktdep = 0.69+1.  
  real :: a12bz,   avgdbz, t69, pmid, mumol, ravgdbz
  real :: vu, kt, df, rhom, tdep
  
  real :: n_o, n_o2, n_o3, n_h2, n_n2, qn2
  real :: cpmult, cvmult,rdmult
  real :: ramo, ramo2, ramo3, ramh2, ramn2, ratio, dpc, dpm 
  real :: vuz(km+1), ktz(km+1), dfz(km+1), rhoz(km+1), wcof, cpx(km+1), dpe(km+1)
  integer  :: i, k
  real :: Runiv
      runiv = 8314.
      a12bz = a12 * bz
      avgdbz= avgd * bz *1000.
      ravgdbz= 1./avgdbz
      ramo =1./amo
      ramo2 =.5*ramo 
      ramo3 =ramo/3. 
      ramh2 = 1./amh2o      
      ramn2 =1./amn2 
      
      do i = is, ie
          dpe(1) = pe(i,1)
       do k=km, 1, -1
          dpe(k+1) = dp2(i,k)        
          qn2 = 1. -qo(i,k)-qo2(i,k)-qo3(i,k)-qh2o(i,k)
          mumol = 1./(qo(i,k)*ramo +qo2(i,k)*ramo2+qo3(i,k)*ramo3+qh2o(i,k)*ramh2+qn2*ramn2)
	  cpx(k+1) = qo(i,k)*cpo +qo2(i,k)*cpo2+qo3(i,k)*cpo3+qh2o(i,k)*cph2o+qn2*cpn2
	  n_o = qo(i,k)* mumol*ramo
	  n_o2= qo2(i,k)*mumol*ramo2
	  n_o3= qo3(i,k)*mumol*ramo3
	  n_h2= qh2o(i,k)*mumol*ramh2
	  n_n2= qn2*mumol*ramn2
          pmid = .5*(pe(i,k)+pe(i,k+1))
          tdep = t2(i,k)  ** ktdep
	  df = a12bz*t2(i,k)**s121/pmid
	  
          vu =  n_o*muo + n_o2*muo2 + n_o3*muo3 + n_n2*mun2 + n_h2*muh2o
          kt =  n_o*lao + n_o2*lao2 + n_o3*lao3 + n_n2*lan2 + n_h2*lah2o	  
          rhom = tdep*runiv/(mumol *pmid) 
	  
       rhoz(k+1) = 	pmid *mumol/t2(i,k)/runiv      !1e-3 * am * plyr(n)/temp(n) / avgdbz 
       vuz(k+1) = vu * rhom 
       ktz(k+1) = kt * rhom 
       dfz(k+1) = df
       
       enddo
       ratio = 1.              !vuz(2)/vuz(3)
       vuz(1) =vuz(2) *ratio 
       ktz(1) =ktz(2) *ratio 
       dfz(1) =dfz(2) *ratio  
          k=1          
          qn2 = 1. -qo(i,k)-qo2(i,k)-qo3(i,k)-qh2o(i,k)
          mumol = 1./(qo(i,k)*ramo +qo2(i,k)*ramo2+qo3(i,k)*ramo3+qh2o(i,k)*ramh2+qn2*ramn2)
          rhoz(1) = pe(i,1) *mumol/t2(i,1)/runiv 
	  cpx(1) = qo(i,k)*cpo +qo2(i,k)*cpo2+qo3(i,k)*cpo3+qh2o(i,k)*cph2o+qn2*cpn2
	  vumol(i,1) = vuz(1)
	  ktmol(i,1) = ktz(1)/cpx(1)
	  dfmol(i,1) = dfz(1)
	  cp_mu(i,k) =cpx(1)
	  do k =1, km         
	    cpmult  =  .5*(cpx(k)+cpx(k+1))
	    vumol(i,k+1) = .5*(vuz(k)+vuz(k+1))
	    ktmol(i,k+1) = .5*(ktz(k)+ktz(k+1))/cpmult
	    dfmol(i,k+1) = .5*(dfz(k)+dfz(k+1))
	    wcof = .5*(rhoz(k)+rhoz(k+1))
	    dpc= .5*(dpe(k)+dpe(k+1))
	    rhomol(i,k+1) =wcof*wcof*grav(i,k+1)/dpc
	    cp_mu(i,k+1) =cpmult
	  enddo   	  
	  k= 1
	  rhomol(i,k) = rhoz(k)*rhoz(k)*grav(i,k)/dp2(i,k)   
      enddo

! print *, 'get_coef_mdif'
 end  subroutine get_coef_mdif 
 
 subroutine get_coef_turb(is, ie, km, t2, u2, v2, dp2, pm2, pe2, wgrav,  cp_mu, vuedd, ktedd, dfedd, con_adj)
 integer , intent(in) :: is, ie, km
 
 logical, intent(in)  ::  con_adj
 real, intent(in), dimension(is:ie,km) :: t2, u2, v2 
 real, intent(in)     ::  dp2(is:ie, km), pm2(is:ie, km), pe2(is:ie, km+1)
  
 real, intent(in), dimension(is:ie,km+1)      :: wgrav,  cp_mu
 real, intent(out), dimension(is:ie,km+1)     :: vuedd, ktedd, dfedd
 integer :: i, k, ii, kk ,j
 
 real :: t(km), u(km), v(km), pe(km), dp(km), grav(km),  cp(km)
 
 vuedd(:,:) = 1.
 ktedd(:,:) = 1.
 dfedd(:,:) = 1.
! do i=is, ie 
  
! enddo
 end subroutine get_coef_turb
!
!================================
!q, pt, w, u, v, dp2, pe2, dpu, dpv,   
!
 subroutine  get_moldiff(dt, vumol, ktmol, dfmol, rhomol, wgrav,   &
                         q, t, w, u, v, dp2, pe2, dpu, dpv, peu, pev, cap3, &
                         j, je, is, ie, isd, ied, jsd, jed, km, nq,  &
	                 ind_h2o, ind_o2, ind_o3p, ind_o3 )
 implicit none
 			 
 integer, intent(in)   :: j, je		  
 integer, intent(in)   :: is, ie, isd,ied, jsd,jed, km, nq
 integer, intent(in)   ::  ind_h2o, ind_o2, ind_o3p, ind_o3 
 
 real, intent(inout)::  q(isd :ied,jsd:jed, km, nq) 
 real, intent(inout)::  t(isd :ied,jsd:jed, km)
 real, intent(in)   ::  cap3(isd :ied,jsd:jed, km)
 real, intent(inout)::  w(isd :ied,jsd:jed, km) 
 real, intent(inout)::  u(isd :ied,jsd:jed+1, km)  
 real, intent(inout)::  v(isd :ied+1,jsd:jed, km)   
 real, intent(in)   ::  dp2(is:ie, km), pe2(is:ie, km+1)
 real, intent(in)   ::  dpu(is:ie, km),dpv(is:ie+1, km)
 real, intent(in)   ::  peu(is:ie, km+1),pev(is:ie+1, km+1) 
 real, intent(in)   ::  dt 

!   real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   !< u-wind (m/s)-DDEBUG=ON
!   real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   !< v-wind (m/s)
 real, intent(inout), dimension(is :ie, km+1) ::  vumol, ktmol, dfmol, rhomol, wgrav
 
 real, dimension(is :ie, km+1) :: vuedd, ktedd, dfedd
 real, dimension(is :ie, km+1) :: cp_mu 
!
!
! 
 real :: q2(is : ie, km), t2(is : ie, km), ut(is:ie, km), vt(is:ie, km)
 real :: pm2(is : ie, km), cap2(is : ie, km) 
 real, dimension(is:ie,km) :: qo, qo2, qh2o, qo3 
 real, dimension(is:ie+1,km) :: v2  
 integer :: itr, i, k, ii, kk, ju
 
 
 if (j == je+1) then 
    ju = je
 else
    ju = j
 endif    
 
 qo(is : ie, 1:km)  = q(is:ie, ju, 1:km, ind_o3p) 
 qo2(is : ie, 1:km) = q(is:ie, ju, 1:km, ind_o2) 
 qh2o(is : ie, 1:km)= q(is:ie, ju, 1:km, ind_h2o) 
 qo3(is : ie, 1:km) = q(is:ie, ju, 1:km, ind_o3) 
 

 wgrav(is:ie, 1:km+1) = 9.8065
 
 do ii = is, ie
 do kk = 1,km
   t2(ii, kk) = t(ii, ju, kk) 
   ut(ii, kk) = u(ii, j, kk)
   vt(ii, kk) = .5*(v(ii, ju, kk) +  v(ii+1, ju, kk))
   cap2(ii,kk) =cap3(ii, ju, kk)
   pm2(ii, kk) = .5*(pe2(ii,kk+1)+pe2(ii,kk))
   
!   (pe2(ii,kk+1)-pe2(ii,kk))/log(pe2(ii,kk+1)/pe2(ii,kk))
 enddo
!     call fv3wam_dry_adjust(km, t2(ii,:), cap2(ii,:), pm2(ii,:), pe2(ii,:), dp2(ii,:))   
 enddo  
 
 call get_coef_mdif(is, ie, km, t2, qo, qo2, qo3, qh2o, dp2, pe2, wgrav,  vumol, ktmol, dfmol, rhomol, cp_mu)
 call get_coef_turb(is, ie, km, t2, ut, vt,  dp2, pm2, pe2, wgrav, cp_mu,  vuedd, ktedd, dfedd, .true.) 
 
! if (is_master() .and. j == je+1) then
!      print *, ' cf_mdif1 ', maxval(wgrav), maxval(vumol), minval(vumol)
!      print *, ' cf_mdif2 ', maxval(ktmol), minval(ktmol)
!      print *, ' cf_mdif3 ', maxval(dfmol), minval(dfmol)   
!      print *, ' cf_mdif4 ', maxval(rhomol), minval(rhomol)                 
! endif

    if (j == je+1) then
    
    t2(is:ie, 1:km) = u(is:ie, j, 1:km) 
     call get_molviscosity(dt, vumol, rhomol, dpu(is:ie,:), peU(is:ie,:), wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Udis') 
    u(is:ie, j, 1:km) = t2(is:ie, 1:km)     
    
    else 
!
! Update all variables j <=je
!    
    call get_molviscosity(dt, ktmol, rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Tdis') 

 do kk = 1,km
 do ii = is, ie
  t(ii, j, kk) = t2(ii, kk) 
!  if (t2(ii,kk) < 90. )  print *, 'fv_molvisZ Tkin', t2(ii,kk), kk
 enddo
 enddo 
!w    
 t2(is:ie, 1:km) = w(is:ie, j, 1:km)
 call get_molviscosity(dt, vumol, rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Wdis') 
 w(is:ie, j, 1:km) = t2(is:ie, 1:km)
 
! u
 t2(is:ie, 1:km) = u(is:ie, j, 1:km) 
 call get_molviscosity(dt, vumol, rhomol, dpu(is:ie,:), peU(is:ie,:), wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Udis') 
 u(is:ie, j, 1:km) = t2(is:ie, 1:km) 
! 
! v
 v2(is:ie+1, 1:km) = v(is:ie+1, j, 1:km) 
! call get_molviscosity_ie1(dt, vumol, rhomol, dpv, peV, wgrav, is, ie, isd, ied, jsd, jed, km, v2, 'Vdis') 
 v(is:ie+1, j, 1:km) = v2(is:ie+1, 1:km) 
 
 call get_molviscosity(dt, dfmol,rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, qo, 'Opdis')
 q(is:ie, j, 1:km, ind_o3p)   = qo(is : ie, 1:km)
 
 call get_molviscosity(dt, dfmol, rhomol,dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, qo2,'O2dis')
 q(is:ie, j, 1:km, ind_o2)   = qo2(is : ie, 1:km) 
 
 call get_molviscosity(dt, dfmol,rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, qo3, 'O3dis')
 q(is:ie, j, 1:km, ind_o3)   = qo3(is : ie, 1:km) 
  
 call get_molviscosity(dt, dfmol, rhomol,dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, qh2o, 'H2dis')
 q(is:ie, j, 1:km, ind_h2o)   = qh2o(is : ie, 1:km)
! print *, 'get_moldiff' 
 
 endif    ! j=je+1
  
 end subroutine get_moldiff
!=============================
!
! 
 subroutine  get_molviscosity(dt, vum, rhomol, dp2, pe, wgrav, is, ie, isd, ied, jsd, jed, km, u,strdis)
 	
  implicit none
  integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km
  real, intent(in)      ::  dt
  real, intent(inout)   ::  u(is :ie, km) 
  real, intent(in)      ::  vum(is:ie, km+1), dp2(is:ie, km), wgrav(is:ie, km+1), rhomol(is:ie, km+1)
  real, intent(in)      ::  pe(is:ie, km+1)
  character(len=*),  intent(in)      ::   strdis  
!
! compute coef for dU/dt=1/rho*d [(rho*K) dU/dz] /dz
!
!   rho*dz = -dp/g   dU/dt = g d[ K rho^2 *g (du/dp)]dp,  rho = Pe/(RT)

  integer, parameter :: dir_fac = 1
  
  integer  :: i, k
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: dta, rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km)
  logical :: mvis_debug
  character(len=128) :: file_dis 
!  
!  
   mvis_debug = .true.
   dta = abs(dt)
   rdt = 1./dta
   
   if((mvis_debug) .and. (is_master())) then
!      print *, 'isd, ied, jsd, jed ', isd, ied, jsd, jed
!      print *, ' is, ie, km ',    is, ie, km 
     file_dis = trim(strdis//'_dedug.form') 
     open(unit=77, file=file_dis, status='unknown', form='formatted')
     write(77,*) is
     write(77,*) ie
     write(77,*) isd, ied, jsd, jed
     write(77,*) km
     write(77,*) vum
     write(77,*) rhomol
     write(77,*) dp2
     write(77,*) pe
     write(77,*) wgrav
     write(77,*) U    
   endif
  do i= is, ie
    ak(1) = 0.
    ck(km) =0.
    dpf = dp2(i,1)
    dpc = pe(i,1)
    gpf = .5*(wgrav(i,1)+wgrav(i,2))
    
    ck(1) = dta*vum(i,2)*gpf*wgrav(i,2)/dpf*rhomol(i,2)
    bk(1) = ck(1) + 1.
    
    dk(1) = u(i,1)
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       gpf = .5*(wgrav(i,k)+wgrav(i,k+1))*vum(i,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*rhomol(i,k) 
       bk(k) = ak(k) + ck(k) + 1.
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       bk(k) = 1. +ak(k) 
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
!         if(mvis_debug) call mpp_error(FATAL,"ERROR MDIF non positive value of plyr") 
   if((mvis_debug) .and. (is_master())) then
    write(77,*) ak
    write(77,*) bk  
    write(77,*) ck
    write(77,*) dk
    write(77,*) ed
    write(77,*) fd    
    write(77,*) U   
    write(77,*) dt     
    close(77)   
!    call mpp_error(FATAL,"STOP in get_molviscosity -DEBUG")                  
   endif  
 end subroutine  get_molviscosity
!=============================================================
!
! 
  subroutine get_molviscosity_ie1(dt, vum, rhomol, dp2, pe, wgrav, is, ie, isd, ied, jsd, jed, km, u, strdis) 
 implicit none
 integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km
 real, intent(in)      ::  dt
 real, intent(inout)   ::  u(is : ie+1, km) 
 real, intent(in)      ::  pe(is:ie+1, km+1), dp2(is:ie+1, km) 
 real, intent(in)      ::  vum(is:ie, km+1),  wgrav(is:ie, km+1), rhomol(is:ie, km+1)
 
 character(len=*),  intent(in)      ::   strdis  

  integer  :: i, k, ii
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km), dta
  logical :: mvis_debug
  integer, parameter :: dir_fac = 1
  character(len=128) :: file_dis 
!  
!  
   mvis_debug = .true.
   dta = abs(dt)
   rdt = 1./dta   
   if((mvis_debug) .and. (is_master())) then

     file_dis = trim(strdis//'_dedug.form') 
     open(unit=77, file=file_dis, status='unknown', form='formatted')
     write(77,*) is
     write(77,*) ie
     write(77,*) isd, ied, jsd, jed
     write(77,*) km
     write(77,*) vum
     write(77,*) rhomol
     write(77,*) dp2(is:ie, 1:km)  
     write(77,*) pe(is:ie, 1:km+1) 
     write(77,*) wgrav
     write(77,*) U(is:ie, 1:km)    
   endif
  do i= is, ie
    ak(1) = 0.
    ck(km) =0.
    dpf = dp2(i,1)
    dpc = pe(i,1)
    gpf = .5*(wgrav(i,1)+wgrav(i,2))
    
    ck(1) = dta*vum(i,2)*gpf*wgrav(i,2)/dpf*rhomol(i,2)
    bk(1) = ck(1) + 1.
    dk(1) = u(i,1)* dta
    
   
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)*dta
       gpf = .5*(wgrav(i,k)+wgrav(i,k+1))*vum(i,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*rhomol(i,k) 
       bk(k) = ak(k) + ck(k) + 1.
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)*dta
       bk(k) = 1. +ak(k)
       
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
          u(i,1) = fd(1)          !    k=1, top-layer
! from top to ground	  
       do k=2,km
          u(i,k) = ed(k)*u(i,k-1) + fd(k)
       enddo
!       
    else               !dirfac - 1
!    
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
!   dt, vum, rhomol, dp2, pe, wgrav, is, ie, isd, ied, jsd, jed, km, u  
!         if(mvis_debug) call mpp_error(FATAL,"ERROR MDIF non positive value of plyr") 
   if((mvis_debug) .and. (is_master())) then
    write(77,*) ak
    write(77,*) bk  
    write(77,*) ck
    write(77,*) dk
    write(77,*) ed
    write(77,*) fd    
    write(77,*) U(is:ie,1:km)   
    write(77,*) dt     
    close(77)   
!    call mpp_error(FATAL,"STOP in get_molviscosity -DEBUG")                  
   endif  
!
! extra-point for V-winds
! 
   ii = ie  
   i = ie+1
   
    ak(1) = 0.
    ck(km) =0.
      
    dpf = dp2(i,1)
    dpc = pe(i,1)
    gpf = .5*(wgrav(ii,1)+wgrav(ii,2))
    
    ck(1) = dta*vum(ii,2)*gpf*wgrav(ii,2)/dpf*rhomol(ii,2)
    bk(1) = ck(1) + 1.
    dk(1) = u(i,1)* dta
    
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)*dta
       gpf = .5*(wgrav(ii,k)+wgrav(ii,k+1))*vum(ii,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*rhomol(ii,k) 
       bk(k) = ak(k) + ck(k) + 1.
    enddo 
       
    k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)* dta
       bk(k) = 1. +ak(k) 
  if (dir_fac == 1 ) then  
       ed(k) =  ak(k)/bk(k)
       fd(k)=   dk(k)/bk(k)         
    do k = km-1, 1, -1
          rinv = 1./(bk(k)-ck(k)*ed(k+1))
          ed(k) = ak(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k+1))*rinv
    enddo   
        u(i,1) = fd(1) 
    do k=2,km
          u(i,k) = ed(k)*u(i,k-1) + fd(k)
    enddo 
    
    else          ! directional factorization: from [top=>bot] or [bot=>top] 
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
    endif          
 end subroutine  get_molviscosity_ie1
! ========================================
!    

 subroutine fv3wam_dry_adjust(levs, t, kappa, prsl, prsi, delp)
!     
! vay oct 2023 for conv-instability in the MLT 
! make dry adjustment of tkin(levs) +> convectively stable configuration
! "warm the conv-unstable layers": relax unstable T-re profile => T-stab
!  Apply GW-scheme only for in the conv-stable atmosphere
!  

      implicit none
      integer             , parameter ::    k_trop = 40       ! SMT-domain
      integer             , parameter ::    niter =  15                
      integer, intent(in)  :: levs
      real, intent(in)    ::  kappa(levs), prsl(levs), delp(levs)
      real, intent(in)    ::  prsi(levs+1)      
      real, intent(inout) ::  t(levs)  
! local
!
!
      real                      ::    gammad, dtdp , zeps, zepsdp, zgamma, rdenom    
      logical                   ::     dodad, ilconv 
      integer                   ::   k, jiter
      real, dimension(levs)     ::   c1dad, c2dad, c3dad, c4dad
      real                      :: dtm, dtins
! from exobase
! to k_trop define the "instability" in the column
!        

       zeps = 2.0e-5     ! set convergence criteria     
       dtdp = (t(levs-1)-t(levs))/(prsl(levs-1)-prsl(levs))
       gammad = .25*(kappa(levs-1)+kappa(levs))/(t(levs-1)+t(levs))/prsi(levs)  
       dodad = (dtdp + zeps) .gt. gammad
!  do i = 1, ncol
!      cappa = 0.5_r8*(cappav(i,2) + cappav(i,1))
!      gammad = cappa*0.5_r8*(t(i,2) + t(i,1))/pint(i,2)
!      dtdp = (t(i,2) - t(i,1))/(pmid(i,2) - pmid(i,1))
!      dodad(i) = (dtdp + zeps) .gt. gammad
!   end do
       
   do k=2, k_trop   
         gammad = .25*(kappa(k-1)+kappa(k))/(t(k-1)+t(k))/prsi(k)
         dtdp = (t(k-1) - t(k))/(prsl(k-1) - prsl(k))
         dodad = dodad .or. (dtdp + zeps).gt.gammad
   end do
   
!         do k = 1, nlvdry
!            c1dad(k) = cappa*0.5_r8*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
!            c2dad(k) = (1._r8 - c1dad(k))/(1._r8 + c1dad(k))
!            rdenom = 1._r8/(pdel(i,k)*c2dad(k) + pdel(i,k+1))
!            c3dad(k) = rdenom*pdel(i,k)
!            c4dad(k) = rdenom*pdel(i,k+1)
!         end do
   
      if( .not. dodad) return
!
! dodat=.T. => make dry adjustment can be "multiple" layers in FV3       
         zeps = 2.0e-5
         do k=1, k_trop
            c1dad(k) = 0.25*(kappa(k+1)+kappa(k))*(prsl(k+1)-prsl(k))/prsi(k)
            c2dad(k) = (1. - c1dad(k))/(1. + c1dad(k))
            rdenom = 1./(delp(k)*c2dad(k) + delp(k-1))
            c3dad(k) = rdenom*delp(k)
            c4dad(k) = rdenom*delp(k+1)
         end do
	 
50       do jiter=1,niter
            ilconv = .true.
            do k=2,k_trop
               zepsdp = zeps*(prsl(k+1) - prsl(k))
               zgamma = c1dad(k)*(t(k+1) + t(k))
	       dtins = zgamma + zepsdp
	       dtm = t(k-1)-t(k)
               if ( dtm >= dtins) then
                  ilconv = .false.
                  t(k) = t(k)*c3dad(k) + t(k+1)*c4dad(k)
                  t(k+1) = c2dad(k)*t(k)
               endif
            end do
            if (ilconv) go to 80           ! convergence => next longitude
         end do
!
! Crude (double) convergence criterion if no convergence in niter iterations
!
          zeps = zeps + zeps
	 
         if (zeps > 1.e-4) then
!            write(6,*)'DADADJ: No convergence in dry adiabatic adjustment'
           else
!            write(6,810) zeps
            go to 50
         endif
      
80    continue
!
!
810   format(//,'FV3WAM DADADJ: Convergence doubled to EPS=',E9.4)
!   COL: do i = 1, ncol
!      if (dodad(i)) then
!         zeps = 2.0e-5_r8
!         do k = 1, nlvdry
!            c1dad(k) = cappa*0.5_r8*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
!            c2dad(k) = (1._r8 - c1dad(k))/(1._r8 + c1dad(k))
!            rdenom = 1._r8/(pdel(i,k)*c2dad(k) + pdel(i,k+1))
!            c3dad(k) = rdenom*pdel(i,k)
!            c4dad(k) = rdenom*pdel(i,k+1)
!         end do
!50       continue
!
!         do jiter = 1, niter
!            ilconv = .true.

!!            do k = 1, nlvdry
!               zepsdp = zeps*(pmid(i,k+1) - pmid(i,k))
!              zgamma = c1dad(k)*(t(i,k) + t(i,k+1))
!
!               if ((t(i,k+1)-t(i,k)) >= (zgamma+zepsdp)) then
!                  ilconv = .false.
!                  t(i,k+1) = t(i,k)*c3dad(k) + t(i,k+1)*c4dad(k)
!                  t(i,k) = c2dad(k)*t(i,k+1)
!               end if
!
!            end do
!            if (ilconv) cycle COL ! convergence => next longitude
!         end do
	     
     return
end subroutine fv3wam_dry_adjust    
 
end module eul_eddy_moldis_mod  
