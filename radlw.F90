
program radlw
!----------------------------------------------------------------------- 
! 
! Purpose: Longwave radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,      only: r8 => shr_kind_r8
use parrrtm,           only: nbndlw, ngptlw
use rrtmg_lw_init,     only: rrtmg_lw_ini
use rrtmg_lw_rad,      only: rrtmg_lw
use radconstants,      only: nlwbands, pcols, ncol, lchnk, pver, pverp, rrtmg_levs, ozone_band
use physconst,         only: cpair
use getdatam,          only: getdata
use planck,            only: integrated_planck

implicit none

   integer, parameter :: ntoplw = 1    ! top level to solve for longwave cooling

!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

!
! Input arguments which are only passed to other routines
!

   real(r8):: pmid(pcols,pver)     ! Level pressure (Pascals)

   real(r8):: aer_lw_abs (pcols,pver,nbndlw) ! aerosol absorption optics depth (LW)

   real(r8):: cld(pcols,pver)      ! Cloud cover
   real(r8):: tauc_lw(nbndlw,pcols,pver)   ! Cloud longwave optical depth by band

!
! Output arguments
!
   real(r8):: qrl (pcols,pver)     ! Longwave heating rate
   real(r8):: qrlc(pcols,pver)     ! Clearsky longwave heating rate
   real(r8):: flns(pcols)          ! Surface cooling flux
   real(r8):: flnt(pcols)          ! Net outgoing flux
   real(r8):: flut(pcols)          ! Upward flux at top of model
   real(r8):: flnsc(pcols)         ! Clear sky surface cooing
   real(r8):: flntc(pcols)         ! Net clear sky outgoing flux
   real(r8):: flutc(pcols)         ! Upward clear-sky flux at top of model
   real(r8):: flwds(pcols)         ! Down longwave flux at surface
   real(r8):: fldsc(pcols)         ! Down longwave clear flux at surface
   real(r8):: fcnl(pcols,pverp)    ! clear sky net flux at interfaces
   real(r8):: fnl(pcols,pverp)     ! net flux at interfaces

   real(r8), pointer, dimension(:,:,:) :: lu ! longwave spectral flux up
   real(r8), pointer, dimension(:,:,:) :: ld ! longwave spectral flux down
   
!
!---------------------------Local variables-----------------------------
!
   integer :: i, k, kk, nbnd         ! indices

   real(r8) :: ful(pcols,pverp)     ! Total upwards longwave flux
   real(r8) :: fsul(pcols,pverp)    ! Clear sky upwards longwave flux
   real(r8) :: fdl(pcols,pverp)     ! Total downwards longwave flux
   real(r8) :: fsdl(pcols,pverp)    ! Clear sky downwards longwv flux

   integer :: inflglw               ! Flag for cloud parameterization method
   integer :: iceflglw              ! Flag for ice cloud param method
   integer :: liqflglw              ! Flag for liquid cloud param method
   integer :: icld                  ! Flag for cloud overlap method
                                 ! 0=clear, 1=random, 2=maximum/random, 3=maximum

   real(r8) :: tsfc(pcols)          ! surface temperature
   real(r8) :: emis(pcols,nbndlw)   ! surface emissivity

   real(r8) :: taua_lw(pcols,rrtmg_levs-1,nbndlw)     ! aerosol optical depth by band

   real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day

   ! Cloud arrays for McICA 
   integer, parameter :: nsubclw = ngptlw       ! rrtmg_lw g-point (quadrature point) dimension
   integer :: permuteseed                       ! permute seed for sub-column generator

   real(r8) :: cicewp(pcols,rrtmg_levs-1)   ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,rrtmg_levs-1)   ! in-cloud cloud liquid water path
   real(r8) :: rei(pcols,rrtmg_levs-1)      ! ice particle effective radius (microns)
   real(r8) :: rel(pcols,rrtmg_levs-1)      ! liquid particle radius (micron)

   real(r8) :: cld_stolw(nsubclw, pcols, rrtmg_levs-1)     ! cloud fraction (mcica)
   real(r8) :: cicewp_stolw(nsubclw, pcols, rrtmg_levs-1)  ! cloud ice water path (mcica)
   real(r8) :: cliqwp_stolw(nsubclw, pcols, rrtmg_levs-1)  ! cloud liquid water path (mcica)
   real(r8) :: rei_stolw(pcols,rrtmg_levs-1)               ! ice particle size (mcica)
   real(r8) :: rel_stolw(pcols,rrtmg_levs-1)               ! liquid particle size (mcica)
   real(r8) :: tauc_stolw(nsubclw, pcols, rrtmg_levs-1)    ! cloud optical depth (mcica - optional)

   ! Includes extra layer above model top
   real(r8) :: uflx(pcols,rrtmg_levs+1)  ! Total upwards longwave flux
   real(r8) :: uflxc(pcols,rrtmg_levs+1) ! Clear sky upwards longwave flux
   real(r8) :: dflx(pcols,rrtmg_levs+1)  ! Total downwards longwave flux
   real(r8) :: dflxc(pcols,rrtmg_levs+1) ! Clear sky downwards longwv flux
   real(r8) :: hr(pcols,rrtmg_levs)      ! Longwave heating rate (K/d)
   real(r8) :: hrc(pcols,rrtmg_levs)     ! Clear sky longwave heating rate (K/d)
   real(r8) :: lwuflxs(nbndlw,pcols,pverp+1)  ! Longwave spectral flux up
   real(r8) :: uoz(nbndlw,pcols,pverp+1)  ! Longwave spectral flux up
   real(r8) :: lwdflxs(nbndlw,pcols,pverp+1)  ! Longwave spectral flux down

   real(r8) :: h2ovmr(pcols,rrtmg_levs) 
   real(r8) :: o3vmr(pcols,rrtmg_levs) 
   real(r8) :: co2vmr(pcols,rrtmg_levs) 
   real(r8) :: ch4vmr(pcols,rrtmg_levs) 
   real(r8) :: o2vmr(pcols,rrtmg_levs) 
   real(r8) :: n2ovmr(pcols,rrtmg_levs) 
   real(r8) :: cfc11vmr(pcols,rrtmg_levs) 
   real(r8) :: cfc12vmr(pcols,rrtmg_levs) 
   real(r8) :: cfc22vmr(pcols,rrtmg_levs) 
   real(r8) :: ccl4vmr(pcols,rrtmg_levs) 
   real(r8) :: pmidmb(pcols,rrtmg_levs) 
   real(r8) :: pintmb(pcols,rrtmg_levs+1) 
   real(r8) :: tlay(pcols,rrtmg_levs) 
   real(r8) :: tlev(pcols,rrtmg_levs+1) 
   real(r8) :: fixflux, tes_reported_o3flux
   !-----------------------------------------------------------------------

   allocate(lu(pcols,rrtmg_levs+1,nlwbands))
   allocate(ld(pcols,rrtmg_levs+1,nlwbands))

   call rrtmg_lw_ini

   cld = 0.
   cicewp = 0.
   cliqwp = 0.
   tauc_lw = 0.
   cld_stolw = 0.
   cicewp_stolw = 0.
   cliqwp_stolw = 0.
   tauc_stolw = 0.
   taua_lw = 0.
   
   icld = 2
   inflglw = 0
   iceflglw = 0
   liqflglw = 0

   call getdata(emis,tsfc,tlev,tlay,pmidmb,pintmb,h2ovmr,o3vmr,co2vmr,ch4vmr,o2vmr,n2ovmr,cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,tes_reported_o3flux)

   if (associated(lu)) lu(1:ncol,:,:) = 0.0_r8
   if (associated(ld)) ld(1:ncol,:,:) = 0.0_r8

   call rrtmg_lw(lchnk  ,ncol ,rrtmg_levs    ,icld    ,                 &
        pmidmb  ,pintmb  ,tlay    ,tlev    ,tsfc    ,h2ovmr, &
        o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,cfc11vmr,cfc12vmr, &
        cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
        cld_stolw,tauc_stolw,cicewp_stolw,cliqwp_stolw ,rei, rel, &
        taua_lw, &
        uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, &
        lwuflxs, lwdflxs, uoz)

   !
   !----------------------------------------------------------------------
   ! All longitudes: store history tape quantities
   ! Flux units are in W/m2 on output from rrtmg_lw and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   ! Heating units are in K/d on output from RRTMG and contain output for
   ! extra layer above model top with vertical indexing from bottom to top.
   ! Heating units are converted to J/kg/s below for use in CAM. 

   flwds(:ncol) = dflx (:ncol,1)
   fldsc(:ncol) = dflxc(:ncol,1)
   flns(:ncol)  = uflx (:ncol,1) - dflx (:ncol,1)
   flnsc(:ncol) = uflxc(:ncol,1) - dflxc(:ncol,1)
   flnt(:ncol)  = uflx (:ncol,rrtmg_levs) - dflx (:ncol,rrtmg_levs)
   flntc(:ncol) = uflxc(:ncol,rrtmg_levs) - dflxc(:ncol,rrtmg_levs)
   flut(:ncol)  = uflx (:ncol,rrtmg_levs)
   flutc(:ncol) = uflxc(:ncol,rrtmg_levs)

   !
   ! Reverse vertical indexing here for CAM arrays to go from top to bottom.
   !
   ful = 0._r8
   fdl = 0._r8
   fsul = 0._r8
   fsdl = 0._r8
   ful (:ncol,pverp-rrtmg_levs+1:pverp)= uflx(:ncol,rrtmg_levs:1:-1)
   fdl (:ncol,pverp-rrtmg_levs+1:pverp)= dflx(:ncol,rrtmg_levs:1:-1)
   fsul(:ncol,pverp-rrtmg_levs+1:pverp)=uflxc(:ncol,rrtmg_levs:1:-1)
   fsdl(:ncol,pverp-rrtmg_levs+1:pverp)=dflxc(:ncol,rrtmg_levs:1:-1)

   
   fnl(:ncol,:) = ful(:ncol,:) - fdl(:ncol,:)
   ! mji/ cam excluded this?
   fcnl(:ncol,:) = fsul(:ncol,:) - fsdl(:ncol,:)

   ! Pass longwave heating to CAM arrays and convert from K/d to J/kg/s
   qrl = 0._r8
   qrlc = 0._r8
   qrl (:ncol,pverp-rrtmg_levs+1:pver)=hr (:ncol,rrtmg_levs-1:1:-1)*cpair*dps
   qrlc(:ncol,pverp-rrtmg_levs+1:pver)=hrc(:ncol,rrtmg_levs-1:1:-1)*cpair*dps

   ! Return 0 above solution domain
   if ( ntoplw > 1 )then
      qrl(:ncol,:ntoplw-1) = 0._r8
      qrlc(:ncol,:ntoplw-1) = 0._r8
   end if

   ! Pass spectral fluxes, reverse layering
   ! order=(/3,1,2/) maps the first index of lwuflxs to the third index of lu.
   if (associated(lu)) then
      lu(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(uoz(:,:ncol,rrtmg_levs:1:-1), &
           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
!      lu(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(lwuflxs(:,:ncol,rrtmg_levs:1:-1), &
!           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
   end if
   
   if (associated(ld)) then
      ld(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(lwdflxs(:,:ncol,rrtmg_levs:1:-1), &
           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
   end if

   call integrated_planck(980._r8, 985._r8, tsfc(1), fixflux)

   !print *, 'tes-reported 03flux','ozone lw up top', 'fixed flux (980-985)'
   print *, tes_reported_o3flux, lu(1,2,ozone_band), lu(1,2,ozone_band)-fixflux
   
end program

