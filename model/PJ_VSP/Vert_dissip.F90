 subroutine get_coef_mdif(is, ie, km, t2, qo, qo2, qo3, qh2o, dp2, pe, grav,  &
                          vumol, ktmol, dfmol, rhomol,cp_mu, am_mol, zgeo)	
 use constants_mod,       only:  rdgas
 implicit none
 
 integer, intent(in)   ::  is, ie , km
 real, dimension(is:ie,km) :: qo, qo2, qh2o, qo3  
 real, dimension(is:ie,km) :: dp2, t2
 real, dimension(is:ie,km+1) :: pe, grav,  vumol, ktmol, dfmol, rhomol
 real, dimension(is:ie,km+1) :: cp_mu, zgeo
 real, dimension(is:ie,km  ) :: am_mol
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
  real :: Runiv, RE, Re2, g981, re2g
      runiv = 8314.
      re = 6370.e3; re2 = re*re
      g981 =9.8065; re2g = re2 * g981
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
	  zgeo(i,km+1) = 0.
       do k=km, 1, -1
          dpe(k+1) = dp2(i,k)        
          qn2 = 1. -qo(i,k)-qo2(i,k)-qo3(i,k)-qh2o(i,k)
          mumol = 1./(qo(i,k)*ramo +qo2(i,k)*ramo2+qo3(i,k)*ramo3+qh2o(i,k)*ramh2+qn2*ramn2)
	  am_mol(i,k) = mumol
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
	  zgeo(i,k) = zgeo(i,k+1) + dp2(i,k)/rhoz(k+1)/g981	
	  grav(i,k) = re2g/(zgeo(i,k)+re)/(zgeo(i,k)+re)
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
	  rhomol(i,k) = 0.25*rhoz(k)*rhoz(k)*grav(i,k)/dp2(i,k)   
      enddo

! print *, 'get_coef_mdif'
 end  subroutine get_coef_mdif 
 
 subroutine get_coef_turb(is, ie, km, t2, u2, v2, dp2, pm2, pe2, gvz, cp_mu, am, zm, vuedd, ktedd, dfedd, &
             con_adj, nfixed, nrest)
 integer , intent(in) :: is, ie, km
 integer ::  nfixed, nrest
 
 logical, intent(in)  ::  con_adj
 real, intent(in), dimension(is:ie,km) :: t2, u2, v2 
 real, intent(in)     ::  dp2(is:ie, km), pm2(is:ie, km), pe2(is:ie, km+1)
 real, intent(in)     ::  am(is:ie, km)  
 real, intent(in), dimension(is:ie,km+1)      ::  gvz,  cp_mu, zm
 real, intent(out), dimension(is:ie,km+1)     :: vuedd, ktedd, dfedd
 integer :: i, k, ii, kk ,j
 
 real :: t(km), u(km), v(km), pe(km), dp(km), grav(km),  kap(km)
 
 real, parameter :: lturb = 35., vumin = 0.03, bn2min = 4.e-8, ric = 0.25
 real, parameter :: dw2min = 0.01, dked_min = 1., dked_max= 1.e6
 real, parameter :: lz_smag = 35., cs_smag = 0.53
 real, parameter :: ruvac = 0.25/8314.
 
 real  :: uz, vz, tz, tc, bn2, dzg, gcp, rik, shr2, ritur
 real  :: bn_uns, zmetk, zlturb, w1, kamp, zgrow, dtemp, ht 
 real  :: t1(km), pt(km), tadj(km), kzz_edd(km+1)
 real  :: lturbmax
 integer  :: ins_con(km), num_ins
 vuedd(:,:) = vumin
 ktedd(:,:) = vumin
 dfedd(:,:) = vumin
 
 nfixed =0
 nrest=0
 lturbmax =300.
 do i=is, ie 
     ins_con(1:km) = 0
     num_ins =0 
     kzz_edd(1:km+1) = vumin
  do k=1, km-1
     dzg =  zm(i,k) - zm(i,k+1) 
     uz =  (u2(i,k)     - u2(i,k+1)) 
     vz =  (v2(i,k)     - v2(i,k+1)) 
     tz =  (t2(i,k)     - t2(i,k+1)) /dzg 
     tc = .5*(t2(i,k)   + t2(i,k+1))
     gcp = gvz(i,k)/cp_mu(i,k)
     bn2 = gvz(i,k)*(tz +gcp)/tc
     shr2 = max(uz*uz+vz*vz, dw2min)
     shr2 = shr2/dzg/dzg
     rik = bn2/shr2
     ht = ruvac/tc*gvz(i,k)*am(i,k)
  if (rik < ric)  then
     if(bn2 < bn2min) then
       bn_uns = bn2
       bn2 = bn2min
	if (bn_uns < 0) then 
	 ins_con(k) = 1
	 num_ins = num_ins + 1
	endif 
     endif  
  endif   
    zgrow = exp(zm(i,k)*Ht)   
    lturbmax = dzg*0.25     
    zlturb = min(lturb*zgrow, lturbmax)
    
           kamp = sqrt(shr2)* zlturb*zlturb
	   ritur = bn2/shr2 *dzg *dzg
           w1 = 1./(1. + 5.*ritur)
	   w1 = w1*w1  
	   	   
!          w1 = 1./(1. + 10.*ritur)/(1.+8.*ritur)  !HB-93 
	            
       kzz_edd(k)= min(max(kamp * w1, dked_min), dked_max)  
  enddo
  
      nfixed = nfixed + num_ins  
!
! check con_adjust_t1   nrest
!      
!    do k=2, km
!       dzg =  zm(i,k) - zm(i,k+1) 
!       dtemp =t2(i,k)     - t2(i,k+1)
!       tz =  dtemp /dzg 
!       tc = .5*(t2(i,k)   + t2(i,k+1))  
!       gcp = gvz(i,k)/cp_mu(i,k)   
!       if ((tz +gcp).lt.0.) nrest= nrest+1       
!    enddo
     
    do k=2,km-1
        w1 = .25*( kzz_edd(k-1)+2.*kzz_edd(k)+kzz_edd(k+1) )     
        ktedd(i,k) =w1
        vuedd(i,k)= w1
        dfedd(i,k)= w1
    enddo	 
 enddo         !i-loop
 
 end subroutine get_coef_turb
 
 subroutine con_adjust_t2( is, ie, km, t2, dp2, pe2, cap2)
 
 
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
 
!====================================================
!q, pt, w, u, v, dp2, pe2, dpu, dpv,   
!
! implicit driver for the vertical M-E dissipation
! + convective adjustment "dry"
!====================================================



 subroutine  get_moldiff(dtin, vumol, ktmol, dfmol, rhomol, wgrav,   &
                         q, t, w, u, v, dp2, pe2, dpu, dpv, peu, pev, cap3, &
                         j, je, is, ie, isd, ied, jsd, jed, km, nq,  &
	                 ind_h2o, ind_o2, ind_o3p, ind_o3 )
 implicit none
 logical  :: con_adj			 
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
 real, intent(in)   ::  dtin 

!   real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   !< u-wind (m/s)-DDEBUG=ON
!   real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   !< v-wind (m/s)

 real, intent(inout), dimension(is :ie, km+1) ::  vumol, ktmol, dfmol, rhomol, wgrav
 
 character(len=128) :: file_dis 
 logical                       :: mvis_debug
 real, dimension(is :ie, km+1) :: vuedd, ktedd, dfedd
 real, dimension(is :ie, km+1) :: cp_mu, zgeo
 real, dimension(is :ie, km)   :: am_mol, kion 
 real :: dt_ic, dt 
 real :: Ne_prof(km)
 real :: nef2, hef2, hch_f2, nede, hkde, hpde, bdde, bede,zchap
!
!

 real :: q2(is : ie, km), t2(is : ie, km), ut(is:ie, km), vt(is:ie, km)
 real :: pm2(is : ie, km), cap2(is : ie, km) 
 real, dimension(is:ie,km) :: qo, qo2, qh2o, qo3 
 real, dimension(is:ie+1,km) :: v2
 real, dimension(is:ie,km)   :: v2c  
 integer :: itr, i, k, ii, kk, ju, nfixed, nrest
!
! km  1/km and cm-3
! 
!  dt_ic = 180.
  dt = dtin
  
  con_adj   =.false.
  mvis_debug=.false.
  
 nef2 = 7.10264e5
 hef2 = 300.e3
 hch_f2 = 39.e3
 nede = 1.43e13
 hkde = 95.e3
 hpde = 70.e3
 bdde = 0.65e-3
 bede =0.14e-3
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
  
 enddo  
  
 if (con_adj) call con_adjust_t2( is, ie, km, t2, dp2, pe2, cap2) 
 
 
 call get_coef_mdif(is, ie, km, t2, qo, qo2, qo3, qh2o, dp2, pe2, wgrav,  vumol, ktmol, dfmol, rhomol, &
        cp_mu, am_mol, zgeo)
	
 do ii = is, ie
  do kk = km, 1,-1  
   zchap = (zgeo(ii,kk)-hef2)/hch_f2 
   Ne_prof(kk) = nef2*exp(0.5*(1- zchap -exp(-zchap)))
   kion(ii,kk) = 7.22e-11* (t2(ii, kk) )** 0.37 *Ne_prof(kk) 
!   kion(ii,kk) = 0.
 enddo 
 enddo	
	
 call get_coef_turb(is, ie, km, t2, ut, vt,  dp2, pm2, pe2, wgrav, cp_mu, am_mol, zgeo, vuedd, ktedd, dfedd,&
                   con_adj, nfixed, nrest) 
   if((mvis_debug) .and. (is_master())) then
!      print *, 'isd, ied, jsd, jed ', isd, ied, jsd, jed
!      print *, ' is, ie, km ',    is, ie, km 
     file_dis = trim('zgeo'//'_dedug.form') 
     open(unit=77, file=file_dis, status='unknown', form='formatted')
     write(77,*) is
     write(77,*) ie
     write(77,*) isd, ied, jsd, jed
     write(77,*) km
     write(77,*) vumol
     write(77,*) vuedd  
     write(77,*) ktmol
     write(77,*) kion     
     write(77,*) zgeo
     write(77,*) cp_mu
     write(77,*) am_mol
     write(77,*) wgrav 
     write(77,*) ut
     write(77,*) vt 
     write(77,*) t2  
     write(77,*)  dp2
     write(77,*)  pm2
     write(77,*)  pe2  
     close(77)                              
   endif  		   
  vumol =vumol +  vuedd
  ktmol =ktmol +  ktedd*0.0
  dfmol =dfmol +  dfedd*0.0
  
 if (is_master() .and. j == je) then
!      print *, ' mdif-kion ', maxval(kion), 1./maxval(kion)*86400.
!      print *, ' mdif-zgeo/grav ', maxval(zgeo)*1.e-3, minval(wgrav)
      print *, ' mdif-ktedd ', maxval(ktedd), minval(ktedd)   
      print *, ' mdif-nfixed ', nfixed, nrest                 
 endif

    if (j == je+1) then
    
    ut(is:ie, 1:km) = u(is:ie, j, 1:km) 
     call uvget_molviscosity(dt, kion, vumol, rhomol, dpu(is:ie,:), peU(is:ie,:), wgrav, is, ie, isd, ied, jsd, jed, km, ut, 'Udis') 
    u(is:ie, j, 1:km) = ut(is:ie, 1:km)     
    
    else 
!
! Update all variables j <=je
!    
    call get_molviscosity(dt, ktmol, rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Tdis') 
 do ii = is, ie
   do kk = 1,km
!
! consider to apply the dry-adjustment for t2(ii, kk)
! 
   t(ii, j, kk) = t2(ii, kk) 
!  if (t2(ii,kk) < 90. )  print *, 'fv_molvisZ Tkin', t2(ii,kk), kk
 enddo
 enddo 
!w    
 t2(is:ie, 1:km) = w(is:ie, j, 1:km)
 call uvget_molviscosity(dt, kion, vumol, rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Wdis') 
 w(is:ie, j, 1:km) = t2(is:ie, 1:km)
 
! u
 t2(is:ie, 1:km) = u(is:ie, j, 1:km) 
 call uvget_molviscosity(dt, kion, vumol, rhomol, dpu(is:ie,:), peU(is:ie,:), wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Udis') 
 u(is:ie, j, 1:km) = t2(is:ie, 1:km) 
! 
! v
 v2c(is:ie, 1:km) = v(is:ie, j, 1:km) 
 call uvget_molviscosity(dt, kion, vumol, rhomol, dpv(is:ie,:), peV(is:ie,:), wgrav, &
       is, ie, isd, ied, jsd, jed, km, v2c, 'Vdis') 
 v(is:ie, j, 1:km) = v2c(is:ie, 1:km) 
 v(ie+1, j, 1:km) = .5*(v2c(ie, 1:km) +v(ie+1, j, 1:km))
 
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
 subroutine  uvget_molviscosity(dt, kion, vum, rhomol, dp2, pe, wgrav, is, ie, isd, ied, jsd, jed, km, u,strdis)
 	
  implicit none
  integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km
  real, intent(in)      ::  dt
  real, intent(inout)   ::  u(is :ie, km)
  real, intent(in)      ::   kion(is:ie, km)
  real, intent(in)      ::  vum(is:ie, km+1), dp2(is:ie, km), wgrav(is:ie, km+1), rhomol(is:ie, km+1)
  real, intent(in)      ::  pe(is:ie, km+1)
  character(len=*),  intent(in)      ::   strdis  
  integer, parameter :: dir_fac = 1
  
  integer  :: i, k
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: dta, rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km)
  logical :: mvis_debug
  character(len=128) :: file_dis 
!  
!  
   mvis_debug = .false.
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
    bk(1) = ck(1) + 1. +kion(i,1)*dta
    
    dk(1) = u(i,1)
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       gpf = .5*(wgrav(i,k)+wgrav(i,k+1))*vum(i,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*rhomol(i,k) 
       bk(k) = ak(k) + ck(k) + 1.+kion(i,k)*dta
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       bk(k) = 1. +ak(k) +kion(i,k)*dta
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
 end subroutine  uvget_molviscosity
 
!=============================================================
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
   mvis_debug = .false.
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
   mvis_debug = .false.
   dta = abs(dt)
   rdt = 1./dta   

   
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
    else               !dirfac = - 1
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


