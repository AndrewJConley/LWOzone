module getdatam


use shr_kind_mod, only: r8 => shr_kind_r8
use parrrtm,           only: nbndlw, ngptlw
use radconstants, only: pcols, rrtmg_levs


implicit none

public :: getdata

save

contains

subroutine getdata(emis,tsfc,tlev,tlay,pmidmb,pintmb,h2ovmr,o3vmr,co2vmr,ch4vmr,o2vmr,n2ovmr,cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr)


   real(r8) :: emis(pcols,nbndlw)   ! surface emissivity
   real(r8) :: tsfc(pcols)          ! surface temperature
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
   real(r8) :: alt(rrtmg_levs)
   real(r8) :: o3irk(rrtmg_levs)

   real(r8) :: r1, r2, r3, r4, r5, r6, r7

   character(len=300) :: junkchar,filename
   real(r8) :: lat, lon, sza, lndfrc, PS, TS, TSerr, emiss, emiss990, emiss1080
   real(r8) :: o3flux, o3column, o3tropcolumn, Troppres, cloudtop, cloudeffod
   integer :: iblnk, line

   integer :: i = 1



   tsfc(:) = 5000.
   tlev(1,:) = 5000.  ! interfaces
   tlay(1,:) = 5000.   ! midpoints
   h2ovmr(1,:) = 1.
   o3vmr(1,:) = 1.
   co2vmr(1,:) = 1.
   ch4vmr(1,:) = 1.8e-6
   o2vmr(1,:) = 0.209
   n2ovmr(1,:) = 1.
   cfc11vmr(1,:) = 1e-40
   cfc12vmr(1,:) = 1e-40
   cfc22vmr(1,:) = 1e-40
   ccl4vmr(1,:) = 1e-40

!   read*,filename

    call get_command_argument(i, filename)
    if (len_trim(filename) == 0) stop

    write (*,*) trim(filename)

   open(unit=15, file=filename, status='old', access='sequential', form='formatted', action='read')
!   open(unit=15, file='/glade/u/home/aconley/ozone/data/TES_r10658_seq114_scn3_R_FM.asc', status='old', access='sequential', form='formatted', action='read')

   read(15,'(A)') junkchar

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) lat

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) lon

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) sza

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=',1.eq.1)
   read(junkchar(iblnk+2:),*) lndfrc

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) PS

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) tsfc

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) TSerr

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) emiss
   emis(:,:nbndlw) = emiss

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) emiss990

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) emiss1080

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) o3flux
   print*,'tes-reported o3 flux',o3flux

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) o3column

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) o3tropcolumn

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) Troppres

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) cloudtop

   read(15,'(A)') junkchar
   iblnk = index(junkchar,'=')
   read(junkchar(iblnk+2:),*) cloudeffod

   read(15,'(A)') junkchar

   read(15,'(A)') junkchar

   read(15,*) pintmb(1,1),r1,tlev(1,1),r2,r3,r4,r5,r6

   do line = 1, rrtmg_levs
     read(15,*) pmidmb(1,line),alt(line),tlay(1,line),h2ovmr(1,line),co2vmr(1,line),n2ovmr(1,line),o3vmr(1,line),o3irk(line)
   end do

   read(15,*) pintmb(1,rrtmg_levs+1),r1,tlev(1,rrtmg_levs+1),r2,r3,r4,r5,r6

   !print*,'pmid',pmidmb
   
   pintmb(1,2:rrtmg_levs) = exp ( (log(pmidmb(1,1:rrtmg_levs-1)) + log(pmidmb(1,2:rrtmg_levs  )) ) / 2._r8 )

   tlev  (1,2:rrtmg_levs) = (tlay  (1,1:rrtmg_levs-1) + tlay  (1,2:rrtmg_levs  ) ) / 2._r8

   !print*,'pintmb',pintmb

   !print*,'tlev',tlev
   !print*,'tlay',tlay



end subroutine getdata

end module
