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

!>@brief fv_deepatm.F90 UFSWAM of
!>       viscosity (mur), conductivity (lam), and diffusivity (d12)
!>       and their efolding time effectiveness
!>@authors V. Yudin and K. Viner (NOAA-SWO/CUA & NOOA-EMC/CPI)
!
! the Deep Atmosphere Formulation (DAF) and transforms for applications of the
!     Shallow Atmosphere (SA) FV3 algorithms, specifically for UFSWAM
!
!   variable gravity, 
!   radial height, 
! & daf_a, daf_du, daf_dv variables to transform
!                  DA-Euler and DA-VL equations in the SA forms
!

module fv_deepatm_mod
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
!
      use constants_mod,      only: grav, rdgas, cp_air,  radius, omega
      use fv_arrays_mod,      only: fv_grid_type, fv_flags_type, fv_atmos_type, fv_grid_bounds_type      
      use tp_core_mod,        only: copy_corners
      use fv_mp_mod,          only: is_master
      implicit none  
      public :: deepatm_init, deepatm_update, deepatm_step_betad
      public :: volume_deep_vs_shal, ufswam_swa_dens
      contains
       
      subroutine deepatm_init(npz, npx, npy, is, ie, js, je, isd, ied, jsd, jed, gridstruct, bd, &
                              phis, delz, grav_int, grav_mid, rad_int, &
                              rdg, daf_a, daf_fcos, ua, va, flagstruct )  
      implicit none
      integer, 	intent(in) :: npz, npx, npy, is, ie, js, je, isd, ied, jsd, jed 
      type(fv_grid_bounds_type), intent(in) :: bd		         
      type(fv_grid_type),  intent(inout), target :: gridstruct
      type(fv_flags_type), intent(inout) :: flagstruct       

      real    :: r_radius2, r_arad2, arad_zmet, grav_dafcor, grav_rad2
      real    :: rgrav
      real    :: ua_int, va_int, dwdt_mid, dwdt_int , wcv, wcu, dt2_rad, wom2     
       
      logical :: flag_geff
      logical :: deep_atmd
           
      real, intent(inout) ::    phis(bd%isd:bd%ied, bd%jsd:bd%jed)  
      real, intent(in)    ::    delz(bd%is:,bd%js:,1:)   
      real, intent(in)    ::    ua(bd%isd:bd%ied,bd%jsd:bd%jed, npz) 
      real, intent(in)    ::    va(bd%isd:bd%ied,bd%jsd:bd%jed, npz)             
                    
      real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   :: grav_mid
      real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1) :: grav_int      
      real,intent(inout),   dimension(bd%is:bd%ie,bd%js:bd%je,npz)      :: rdg        !rdgas/gravity 
      real, dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1)                :: rad_int       
      
      real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed, npz)   :: daf_a                    
      real, intent(out), dimension(bd%isd:bd%ied, bd%jsd:bd%jed)         :: daf_fcos

      integer :: i, j, k 
      
      
      r_arad2  = 0.5/radius 
      rgrav = 1./grav
      
      grav_rad2=grav*radius*radius           
         
      grav_int(:,:,:) = grav
      grav_mid(:,:,:) = grav
 
      rdg(:,:,:) = -rdgas*rgrav  
      
      flag_geff   = flagstruct%ufswam_geff
      
      deep_atmd =  flagstruct%deep_atmd 
!
! just in case copy_corners of phis to define phis_u/phis_v 
!   
           call copy_corners(phis, npx, npy, 1, gridstruct%nested, bd, &
                            gridstruct%sw_corner, gridstruct%se_corner, &
                            gridstruct%nw_corner, gridstruct%ne_corner)  
!			    
      if(.not.deep_atmd.and.flagstruct%var_grav) then
        do j=js,je
          do i=is,ie
            rad_int(i,j,npz+1) = radius + phis(i,j)*rgrav
            grav_int(i,j,npz+1) = grav_rad2/rad_int(i,j,npz+1)/rad_int(i,j,npz+1)
            do k=npz,1,-1
              rad_int(i,j,k) = rad_int(i,j,k+1) - delz(i,j,k)
              grav_int(i,j,k) =  grav_rad2/rad_int(i,j,k)/rad_int(i,j,k)
              grav_mid(i,j,k) = .5*(grav_int(i,j,k+1)+grav_int(i,j,k))
              rdg(i,j,k) = -rdgas/grav_mid(i,j,k)
	      daf_a(i,j,k) = 1.0
            enddo
          enddo
        enddo
      endif     
       
      if(deep_atmd) then
!  
! ua/va -defined  or call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd, gridstruct%bounded_domain)
! in ROT3/do_f3d: gridstruct%w00(i,j) = 2.d0*omega*cos(agrid(i,j,2)) <........> daf_fcos(i,j)
!
	 ua_int =0.0
	 va_int =0.0
 	
        do j=js,je	  
          do i=is,ie
	   daf_fcos(i,j) =2.*omega*gridstruct%cosa_s(i,j)	  			
	 enddo 
	enddo
		
        do j=js,je	  
          do i=is,ie	
!
! cell-centered A-grid w-equations dP/dr + dens*grav_var_int  =0
!	   	       
	    arad_zmet =  radius + phis(i,j)*rgrav	    
            rad_int(i,j,npz+1)  = arad_zmet
            grav_int(i,j,npz+1) = grav_rad2/arad_zmet/arad_zmet
	    		    
           do k=npz,2,-1
	      arad_zmet = rad_int(i,j,k+1) - delz(i,j,k)
              rad_int(i,j,k) = arad_zmet
!	      
! interface correction-factor for grav_var_int with [ua=> ua_int va => va_int]
            if(flag_geff) then	
	      ua_int = .5*(ua(i,j,k)+ ua(i,j,k-1))
	      va_int = .5*(va(i,j,k)+ va(i,j,k-1))	           
              grav_dafcor= (ua_int*ua_int+va_int*va_int)/arad_zmet +ua_int*daf_fcos(i,j)  
	    else
	      grav_dafcor=0.
	    endif   	     
              grav_int(i,j,k) = grav_rad2/arad_zmet/arad_zmet - grav_dafcor	      
              grav_mid(i,j,k) = .5*(grav_int(i,j,k+1)+grav_int(i,j,k))
              rdg(i,j,k) = -rdgas/grav_mid(i,j,k)
! (Re+Z)/Re
!       [interface(k+1)+interface(k)]/2R = 1. + Z/Re
!	      
	      daf_a(i,j,k) = (arad_zmet+rad_int(i,j,k+1))*r_arad2 	      
            enddo 
! model-top	    
	      k= 1
	      arad_zmet = rad_int(i,j,k+1) - delz(i,j,k)	
	      rad_int(i,j,k) = arad_zmet
           if(flag_geff) then	 
	      ua_int = ua(i,j,k)
	      va_int = va(i,j,k)	      
              grav_dafcor= (ua_int*ua_int+va_int*va_int)/arad_zmet +ua_int*daf_fcos(i,j)
	    else
	      grav_dafcor= 0.
	    endif         
              grav_int(i,j,k) = grav_rad2/arad_zmet/arad_zmet - grav_dafcor	      
              grav_mid(i,j,k) = .5*(grav_int(i,j,k+1)+grav_int(i,j,k))
              rdg(i,j,k) = -rdgas/grav_int(i,j,k)	      
	      daf_a(i,j,k) = (arad_zmet+rad_int(i,j,k+1))*r_arad2  ! [(R+z)/R] at mid-points sum(delz) = Zint              	    
          enddo
        enddo 

       endif
       return		 
!==================================================================            
     end subroutine deepatm_init 
     
     subroutine deepatm_update(npz, is, ie, js, je, isd, ied, jsd, jed, bd,   &
                              phis, delz, grav_int, grav_mid, rad_int, &
			      rdg, daf_a, ua, va, daf_fcos, flagstruct ) 
      implicit none 		      
      integer, 	intent(in) :: npz, is, ie, js, je, isd, ied, jsd, jed 
      type(fv_grid_bounds_type), intent(IN) :: bd				      
!      type(fv_grid_type),  intent(inout), target :: gridstruct
      type(fv_flags_type), intent(inout) :: flagstruct       
      
  
      real    :: r_arad2, arad_zmet, grav_dafcor, grav_rad2
      real    :: rgrav, arad_zmetU, arad_zmetV, delz_u, delz_v, aradU, aradV
      real    :: ua_int, va_int, dwdt_mid, dwdt_int , wcv, wcu, dt2_rad, wom2 
      real    :: rho_dzpt   
       
      logical :: flag_geff
      logical :: deep_atmd
      
      real, intent(in)    ::    delz(bd%is:,bd%js:,1:)        
      real, intent(in)    ::    phis(bd%isd:bd%ied,bd%jsd:bd%jed)  
      real, intent(in)    ::    daf_fcos(bd%isd:bd%ied,bd%jsd:bd%jed) 

      real, intent(in)    ::    ua(bd%isd:bd%ied,bd%jsd:bd%jed, npz) 
      real, intent(in)    ::    va(bd%isd:bd%ied,bd%jsd:bd%jed, npz) 
                                         
      real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz)   :: grav_mid
      real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1) :: grav_int      
      real,intent(inout),   dimension(bd%is:bd%ie,bd%js:bd%je,npz)      :: rdg        !rdgas/gravity 
      real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1) :: rad_int       

      real, intent(inout), dimension(bd%isd:bd%ied,bd%jsd:bd%jed, npz)  :: daf_a     
 
      integer :: i, j, k 
      r_arad2  =0.5/radius 
      rgrav = 1./grav    
      grav_rad2=grav*radius**2 
     
      deep_atmd =  flagstruct%deep_atmd 
      flag_geff  =  flagstruct%ufswam_geff 
      
!!!!!!$OMP parallel do default(none) shared(is,ie,js,je, npz, ps,delp,grav_var, grav_var_int, phis,arad_zmet, &
!!!!!$OMP           rad_int, rdg, grav_rad2, delz)      
        
      if(.not.deep_atmd.and.flagstruct%var_grav) then
        do j=js,je
          do i=is,ie
            rad_int(i,j,npz+1) = radius + phis(i,j)*rgrav
            grav_int(i,j,npz+1) = grav_rad2/rad_int(i,j,npz+1)/rad_int(i,j,npz+1)
            do k=npz,1,-1
              rad_int(i,j,k) = rad_int(i,j,k+1) - delz(i,j,k)
              grav_int(i,j,k) =  grav_rad2/rad_int(i,j,k)/rad_int(i,j,k)
              grav_mid(i,j,k) = .5*(grav_int(i,j,k+1)+grav_int(i,j,k))
              rdg(i,j,k) = -rdgas/grav_mid(i,j,k)
	      daf_a(i,j,k) = 1.
            enddo
          enddo
        enddo
      endif     
       
      if(deep_atmd) then
!  
! ua/va -defined  call c2l_ord2(u, v, ua, va, gridstruct, npz, gridstruct%grid_type, bd, gridstruct%bounded_domain)
! 
	 ua_int =0.0
	 va_int =0.0
!!!!!$OMP parallel do default(none) shared(is,ie,js,je, npz, ps,delp,grav_var, grav_var_int, daf_fcos, phis,arad_zmet, &
!!!!!$OMP           rad_int, rdg, grav_rad2, delz, ua_int,va_int, ua, va, grav_dafcor, flag_geff)  
        do j=js,je	  
          do i=is,ie	
!
! cell-centered A-grid w-equations dP/dr + dens*grav_var_int  =0
!	   	       
	    arad_zmet =  radius + phis(i,j)*rgrav  
            rad_int(i,j,npz+1) = arad_zmet
            grav_int(i,j,npz+1) = grav_rad2/arad_zmet/arad_zmet
	    		    
           do k=npz,2,-1
	      arad_zmet = rad_int(i,j,k+1) - delz(i,j,k)
!	      
! interface correction-factor for grav_var_int with [ua=> ua_int va => va_int]
            if(flag_geff) then	
	      ua_int = .5*(ua(i,j,k)+ ua(i,j,k-1))
	      va_int = .5*(va(i,j,k)+ va(i,j,k-1))	           
              grav_dafcor= (ua_int*ua_int+va_int*va_int)/arad_zmet +ua_int*daf_fcos(i,j)  
	    else
	      grav_dafcor=0.
	    endif   	     
              grav_int(i,j,k) = grav_rad2/arad_zmet/arad_zmet - grav_dafcor	      
              grav_mid(i,j,k) = .5*(grav_int(i,j,k+1)+grav_int(i,j,k))
              rdg(i,j,k) = -rdgas/grav_mid(i,j,k)
! (Re+Z)/Re
!       [interface(k+1)+interface(k)]/2R = 1. + Z/Re
!	      
	      daf_a(i,j,k) = (arad_zmet+rad_int(i,j,k+1))*r_arad2 
	      rad_int(i,j,k) = arad_zmet	      
            enddo 
! model-top	    
	      k= 1
	      arad_zmet = rad_int(i,j,k+1) - delz(i,j,k)	
	      rad_int(i,j,k) = arad_zmet
           if(flag_geff) then	 
	      ua_int = ua(i,j,k)
	      va_int = va(i,j,k)	      
              grav_dafcor= (ua_int*ua_int+va_int*va_int)/arad_zmet +ua_int*daf_fcos(i,j)
	    else
	      grav_dafcor= 0.
	    endif         
              grav_int(i,j,k) = grav_rad2/arad_zmet/arad_zmet - grav_dafcor	      
              grav_mid(i,j,k) = .5*(grav_int(i,j,k+1)+grav_int(i,j,k))
              rdg(i,j,k) = -rdgas/grav_mid(i,j,k)
	      
	      daf_a(i,j,k) = (arad_zmet+rad_int(i,j,k+1))*r_arad2  ! [(R+z)/R] at mid-points sum(delz) = Zint               	    
          enddo
        enddo     
     endif	  
!               			      
     end subroutine deepatm_update
!     
     subroutine deepatm_step_betad(npz, is, ie, js, je, isd, ied, jsd, jed, bd, bdt, &
                u, v, w, daf_a, rad_int, gridstruct, flagstruct)
!
! add impacts of the beta-d terms to the momentum u-v-w after VL-remap on EL-layers
! energy conserving "simultaneous time-advance to conserve kinetic energy due to 
!  cancellation of beta-d terms in Ekin-budget 
!        
!  
      implicit none
      type(fv_grid_bounds_type), intent(IN) :: bd
      type(fv_grid_type),  intent(in), target :: gridstruct
      type(fv_flags_type), intent(in) :: flagstruct    
              
      integer, intent(in) :: npz, is, ie, js, je, isd, ied, jsd, jed
      real, intent(in)    :: bdt      

      real, intent(inout)    ::    w(bd%isd:bd%ied,bd%jsd:bd%jed, npz) 
      real, intent(inout)    ::    v(bd%isd:bd%ied+1,bd%jsd:bd%jed, npz)  
      real, intent(inout)    ::    u(bd%isd:bd%ied,bd%jsd:bd%jed+1, npz) 
      
      real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz)    :: daf_a                        
      real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1)  :: rad_int    

                
! locals   
      integer :: i, j, k     
      real :: bdt2, wom2, wcu, ehkin, rad_mid, wombdt, dt2_rad, rdafu
      
       dt2_rad = bdt/radius
       wom2 =2.*omega
       bdt2 =0.5*bdt
       
     if(.not.flagstruct%ufswam_geff) then
      
      do k=1,npz
        do j=js,je
         do i=is,ie  
	  rad_mid = rad_int(i,j,k)+rad_int(i,j,k+1)
	  wcu = u(i,j,k)+u(i,j+1,k)
	  ehkin = wcu*wcu+(v(i,j,k)+v(i+1,j,k))**2
!	  
! 0.25*2 = 0.5*bdt ,ie+1 for V je+1 for U cell-centered values + mid-layer 1./radius
!  
	 w(i,j,k) = w(i,j,k) + bdt2*(ehkin/rad_mid + wom2*gridstruct%cosa_s(i,j)*wcu)
         enddo
        enddo
       enddo	
       	    	
     endif
!      
!  allocate ( cosa_s(isd_2d:ied_2d,jsd_2d:jed_2d) )  ! cell center for w
!  allocate   cosa_u(isd_2d:ied_2d+1,jsd_2d:jed_2d)   I+1 for v_dgrid
!  allocate ( cosa_v(isd_2d:ied_2d,jsd_2d:jed_2d+1) ) J+1 for u_dgrid
!   
!  use daf_dv daf_du to get D-rad_int terma after the complete bdt-time_step
!   
       wombdt  =   wom2*bdt
      do k=1,npz
        do j=js,je+1
         do i=is,ie
          wcu =0.5*(w(i,j,k)+w(i,j-1,k))
	  rdafu =1./(daf_a(i,j,k)+daf_a(i,j-1,k))
          u(i,j,k) = u(i,j,k) -  wcu*(dt2_rad*u(i,j,k)*rdafu+wombdt*gridstruct%cosa_v(i,j))
         enddo
        enddo
!	
        do j=js,je
          do i=is,ie+1
           wcu =(w(i,j,k)+w(i-1,j,k))/(daf_a(i,j,k)+daf_a(i-1,j,k))
           v(i,j,k) = v(i,j,k) - dt2_rad*v(i,j,k)*wcu
         enddo 
	 enddo  
       enddo   
    

     end subroutine deepatm_step_betad
     
     subroutine volume_deep_vs_shal(npz, is, ie, js, je, isd, ied, jsd, jed, bd, &
                rad_int, delz, gridstruct, flagstruct) 
      implicit none
      type(fv_grid_bounds_type), intent(IN) :: bd
      type(fv_grid_type),  intent(in), target :: gridstruct
      type(fv_flags_type), intent(in) :: flagstruct    
              
      integer, intent(in) :: npz, is, ie, js, je, isd, ied, jsd, jed
    

      real, intent(in) ::   delz(bd%is:,bd%js:,1:)                      
      real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed,npz+1)  :: rad_int 
! local
      real, dimension(bd%is:bd%ie,bd%js:bd%je,npz) :: vol_da, vol_sa
      real :: tvol_da, tvol_sa, rad_mid,rda2, vda_cor , dza
      
      integer :: i, j, k   
       tvol_da =0.
       tvol_sa =0.
        do k=1,npz
         do j=js,je
          do i=is,ie  
	  dza = abs(delz(i,j,k))
	  rad_mid =.5*(rad_int(i,j,k+1)+ rad_int(i,j,k))
! area corr	  
	  rda2 = (rad_mid/radius)**2
! radial corr	  
	  vda_cor = 1. + dza/rad_mid +1./3.*(dza/rad_mid)**2.
	  
	   vol_sa(i,j,k) = gridstruct%area(i,j)*dza
	   vol_da(i,j,k) = gridstruct%area(i,j)*rda2*dza*vda_cor 
	   tvol_da = tvol_da +vol_da(i,j,k)
	   tvol_sa = tvol_sa +vol_sa(i,j,k)
          enddo 
	 enddo  
       enddo  	
        if (is_master () ) write(*,*) ' volume_deep_vs_shal ',  tvol_da,  tvol_sa,  tvol_da/tvol_sa
	
!!!! 0:   volume_deep_vs_shal   3.8431319E+17  3.5120391E+17   1.094274
			 
     end subroutine volume_deep_vs_shal
!     
     subroutine ufswam_swa_dens(npz, is, ie, js, je, isd, ied, jsd, jed, bd,   &
                pt, pe, peln, delp, delz, rd_msp, grav_var, wam_dens,wam_dhyd, &
	        wam_zhyd, phis, gridstruct, flagstruct)  
		
      implicit none
      type(fv_grid_bounds_type), intent(IN) :: bd
      type(fv_grid_type),  intent(in), target :: gridstruct
      type(fv_flags_type), intent(in) :: flagstruct    
              
      integer, intent(in) :: npz, is, ie, js, je, isd, ied, jsd, jed
    

      real, intent(in)    ::   delz(bd%is:,bd%js:,1:)  
      real, intent(in)    ::   phis(bd%isd:bd%ied,bd%jsd:bd%jed)
      real, intent(inout) ::  peln(bd%is:bd%ie,npz+1,bd%js:bd%je)           !< ln(pe)                        
      real, intent(inout) ::  pe(bd%is-1:bd%ie+1, npz+1,bd%js-1:bd%je+1)  !< edge pressure (pascal)
      
      real, intent(in)    ::   pt(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      real, intent(in)    ::   delp(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      real, intent(in)    ::   grav_var(bd%isd:bd%ied,bd%jsd:bd%jed,npz) 
      real, intent(in)    ::   rd_msp(bd%isd:bd%ied,bd%jsd:bd%jed,npz)       
      real, intent(out)   ::   wam_dens(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      real, intent(out)   ::   wam_dhyd(bd%isd:bd%ied,bd%jsd:bd%jed,npz)
      real, intent(out)   ::   wam_zhyd(bd%isd:bd%ied,bd%jsd:bd%jed,npz)                            	
! local

      real   ::   wzmhyd, wzkhyd,  wdenshyd 
      real   ::   wam_znh(bd%is:bd%ie,bd%js:bd%je,npz), zhm, zhk  
      integer :: i, j, k        
!!!$OMP parallel do default(none) shared(is,ie,js,je,npz, num_gas, wzmhyd,wzkhyd, wam_dens, wam_dhyd, wzhyd, &
!!!$OMP                           phis, pt, pe, peln,delp, delz, grav_var)      
      
        do j=js,je
          do i=is,ie     
	    wzmhyd  = phis(i,j)/grav
	    zhm =wzmhyd
	   do k=npz, 1, -1
! hydrostatic height 	  wzkhyd > 0 zm_hyd vs zm_nh   g=const  549.  vs nh-delz 602. dz_top= ~50/60 km
!
	     wzkhyd = wzmhyd + rd_msp(i,j,k)*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))/grav_var(i,j,k)
	     
!	     wzkhyd = wzmhyd + rd_msp(i,j,k)*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))/grav

	     wam_zhyd(i,j,k) = .5*(wzkhyd +wzmhyd)
	     wzmhyd =wzkhyd
	     
	     zhk = zhm - delz(i,j,k)
	     wam_znh(i,j,k) = .5*(zhk+zhm)
	     zhm = zhk
	     
! hydrostatic density 	     
	     wdenshyd = -delp(i,j,k)/delz(i,j,k)/grav_var(i,j,k) 
! state-law:  wdens =P/(RT)
	     wam_dens(i,j,k) = .5*(pe(i,k, j)+pe(i,k+1, j))/pt(i,j,k)/rd_msp(i,j,k)
             wam_dhyd(i,j,k) = wdenshyd 
	   enddo 	   
	  enddo 
	enddo  
    if (is_master()) then 
	write(*, *), 'dens_wam ', maxval(wam_dens(is:ie, js:je,1:npz)), minval(wam_dens(is:ie, js:je,1:npz)) 
	write(*, *), 'dens_hyd ', maxval(wam_dhyd(is:ie, js:je,1:npz)), minval(wam_dhyd(is:ie, js:je,1:npz)) 
	write(*, *), 'zm_hyd vs zm_nh ', maxval(wam_zhyd(is:ie, js:je,1:npz)), maxval(wam_znh(is:ie, js:je,1:npz)) 
	
    endif		      	
     end subroutine ufswam_swa_dens		     			                                       
end module fv_deepatm_mod        
