! 2005/11/17 ��������
! 2006/05/09 �ǽ�����
! �̳�ƻ����������ϵ�ʳز�
! ����ʪ���ظ��漼��ǯ
! ��¼ �� (Yumura Tsubasa)
! �쥤�꡼�ƥ��顼�԰���

!####--- �ѥ�᡼�������� ---####
program main
implicit none

integer x,z,t,xmax,zmax,tmax,datanum,i,j,imax,ic,hassant,hassan, &
& jtype,jadd,Ttype,exisped,exiswind,bdet,datanumn,ntype
character filename*50
character filename2*50
character filename3*50
character filename4*50
real*8, dimension(-1:129,-1:129) :: n,n00,dataTemp,phi,phiold,phinew, &
& vix,viz,vex,vez,jx,jz,divj,divE,U,Eped,datan,pedersen,wind,nvx,nvz
real*8, dimension(0:101) :: Tempn,Tempi,HO,n1,nneu,nu,gamma,ki,bi, &
& alpha,beta,bidash,E0
real*8 dt,zbot,ztop,xscale,zscale,zpeak,npeak,chi,dx,dz, &
& mO,mN,g,k,datazmin,datazmax,datadz,height,e,B,omega,Temp0,nneu0, &
& z1,dataz,A,ni,zz,dn,jouran,dnx,dnz,dphix,dphiz,phisumx,phisumz,jrate, &
& dvix,dviz,dvex,dvez,djx,djz,U0,n0,dnvx,dnvz,difn,sumdifn,bheight,bpoint

zpeak=3.0e5            ! �ץ饺�ޥԡ����ι��� [m]
chi=60.0/90.0*3.14/2.0 ! ����ŷĺ�Ѧ� [rad]
jtype=1                ! �������η�, 0:̵��, 1:sin��, 2:sin�ȿ�3, 3:����, 4:ˢ
jadd = 1               ! ��������Ϳ����, 1:̩��, 2:���� 
jrate=0.01             ! �������γ��
Ttype=1                ! ���٤�Ϳ����, 1:��, 2:��, 3:ͼ��
ntype=3                ! ̩�٤�Ϳ����, 1:��, 2:��, 3:�� 
exisped=0              ! �ڥ���������ή��̵ͭ, 0:̵��, 1:ͭ��
exiswind=0             ! ��������̵ͭ,         0:̵��, 1:ͭ��

tmax=20001             ! �����ॹ�ƥå׿�
imax=1000              ! ��«
dt=1.0                 ! ���ֳִ� [s]
zbot=1.0e5             ! �ǲ����ι��� [m]
ztop=6.0e5             ! �Ǿ����ι��� [m]
xmax=100                ! x �γʻҿ�
zmax=100                ! z �γʻҿ�
xscale=5.0e4           ! ���٤�Ĺ���������� [m]
zscale=ztop-zbot       ! ���٤�Ĺ���������� [m]


!E0=0.0
U0=100.0              ! ������® [m/s]
npeak=1.0e12           ! �ץ饺�ޥԡ����ο�̩�� [1/m^3]
n0=0.0
dx=xscale/dble(xmax)   ! ���٤γʻҴֳ� [m]
dz=zscale/dble(zmax)   ! ���٤γʻҴֳ� [m]
mO=16.*1.672e-27       ! ������ 16��+ �μ��� [kg]
mN=28.*1.672e-27       ! ����ʬ�� N2 �μ��� [kg]
g=9.8                  ! ���ϲ�®�� [m/s^2]
k=1.381e-23            ! �ܥ�ĥޥ���� [J/K]
e=1.6e-19              ! �Ų� [C]
B=4.5e-5               ! ���� [T]
omega=e*B/mO           ! ������β�ž���ȿ� [1/s]
Temp0=288.             ! ��ɽ�̲��� [K]
nneu0=1.0e24           ! ��ɽ�̤Ǥ������絤��̩�� [1/m^3]
bpoint=8.5e11

! --���ٹ�¤�Υơ��֥��ɤ߹���--
  datanum=69
  if (Ttype==1) then
    open(20,file='dataTemp-night.d')
  else if (Ttype==2) then
    open(20,file='dataTemp-day.d')
  else
    open(20,file='dataTemp-eve.d')
  endif
  read(20,*) ((dataTemp(i,j),i=1,4),j=1,datanum)
  close(20)
  datazmin=dataTemp(1,1)
  datazmax=datatemp(1,datanum)

! --̩�ٹ�¤�Υơ��֥��ɤ߹���--
  datanumn=100
  if (ntype==1) then
    open(21,file='n-high.d')
  else if (ntype==2) then
    open(21,file='n-mid.d')
  else
    open(21,file='n-low.d')
  endif
  read(21,*) ((datan(i,j),i=1,2),j=1,datanum)
  close(21)
  datazmin=datan(1,1)
  datazmax=datan(1,datanum)


!##### ���ٰ�¸����ѥ�᡼�������� #####
do z=1,zmax+1
  height=zbot+dz*dble(z)

!-- �絤�β��� [K] �η׻� --
    do i=2,datanum
      if (dataTemp(1,i) > height) then
        Tempn(z)=dataTemp(2,i-1) &
&            +(height-dataTemp(1,i-1)) &
&            *(dataTemp(2,i)-dataTemp(2,i-1)) &
&            /(dataTemp(1,i)-dataTemp(1,i-1))
        goto 30
      endif
    enddo

30 continue

! -- ������β��� [K] �η׻� --
    do i=2,datanum
      if (dataTemp(1,i) > height) then
        Tempi(z)=dataTemp(3,i-1) &
            +(height-dataTemp(1,i-1)) &
&            *(dataTemp(3,i)-dataTemp(3,i-1)) &
&            /(dataTemp(1,i)-dataTemp(1,i-1))
        goto 31
      endif
    enddo

31  continue

! -- ��������ϥ��� [m] �η׻� --
  HO(z)=k*Tempi(z)/mO/g ! O+ �Υ�������ϥ���
!  HN=k*Tempn(z)/mN/g ! N2 �Υ�������ϥ���

! -- ������ο�̩�� [1/m^3] �η׻� --
  z1=(dble(z)*dz-zpeak+zbot)/HO(z)
!  n1(z)=npeak*cos(chi)*exp(1.0-z1-exp(-z1)/cos(chi))

    do i=2,datanumn
      if (datan(1,i) > height) then
        n1(z)=1.0e6*(datan(2,i-1) &
&            +(height-datan(1,i-1)) &
&            *(datan(2,i)-datan(2,i-1)) &
&            /(datan(1,i)-datan(1,i-1)))

        goto 32
      endif
    enddo

32 continue

! -- �����絤�ο�̩�� [1/m^3] �η׻� --
  if (z==1) then
    nneu(z)=Temp0/Tempn(z)*nneu0 &
&           *exp(-0.5*mN*g/k*height*(1.0/Tempn(z)+1.0/Temp0))
  else
    nneu(z)=Tempn(z-1)/Tempn(z)*nneu(z-1) &
&           *exp(-0.5*mN*g/k*dz*(1.0/Tempn(z)+1.0/Tempn(z-1)))
  endif

! -- ������������絤��ʿ��ʬ���� [g/mol] �η׻� --
  A=(16.*n1(z)+28.*nneu(z))/(n1(z)+1.0e10+nneu(z))

! -- ������������絤�ξ��ͼ��ȿ� [1/s] �η׻� --
  nu(z)=2.6*(n1(z)+1.0e10+nneu(z)*10.0)/1.0e15/A**0.5

! -- ��ž���ȿ��Ⱦ��ͼ��ȿ�����η׻� --
  ki(z)=e*B/mO/nu(z)

! -- ������β�ư���η׻� --
  bi(z)=e/mO/nu(z)

! -- �� = 1/ki*2 �η׻� --
  alpha(z)=1.0/ki(z)**2

! -- E0
  if (exisped==1) then
    E0(z)=g*B/nu(z)
  else
    E0(z)=0.0d0
  endif

enddo

! -- �� = ��(1/ki**2)/��z �η׻� --
ki(0)=ki(1)
do z=1,zmax
  beta(z)=-(ki(z+1)-ki(z-1))/ki(z)**3/dz
enddo

! -- bi' = ��bi/��z �η׻�
do z=1,zmax
  bidash(z)=(bi(z+1)-bi(z-1))/2.0/dz
enddo

!-- �԰�������ĹΨ [1/s] �η׻�
open(19,file='/home/yumu/data/gamma.d')
do z=2,zmax-1
  gamma(z)=(g/nu(z)+E0(z)/B)/n1(z)*(n1(z+1)-n1(z-1))/2.0
  write(19,*) zbot+dz*dble(z),gamma(z)
enddo
close(19)

!####--- ����ͤ����� ---####
!-- n �ν���� --
do x=1,xmax
  do z=1,zmax

!--�������η��λ��� --
!0:̵��
    if (jtype==0) then
      jouran=0.0
    endif

!1:sin ��
    if (jtype==1) then
      jouran=cos(2.0*3.14*dble(x)/dble(xmax))
    endif

!2:sin ��-�ȿ�3
    if (jtype==2) then
      jouran=cos(3.0*2.0*3.14*dble(x)/dble(xmax))
    endif

!3:����
    if (jtype==3) then
      if (x < 0.5*xmax) then
        jouran=-2.0*dble(x)/dble(xmax)+0.5
      else
        jouran=2.0*dble(x)/dble(xmax)-0.5
      endif
    endif

!4:ˢ
    if (jtype==4) then
      if ((x-xmax/4.0)**2+(z-zmax/3)**2 < 9) then
        jouran=-1.0
        jrate=0.05
      else if ((x-xmax/2.0)**2+(z-zmax/3)**2 < 9) then
        jouran=-1.0
        jrate=0.05
      else if ((x-xmax*3.0/4.0)**2+(z-zmax/3)**2 < 9) then
        jouran=-1.0
        jrate=0.05
      else
        jouran=0.0
      endif
    endif


!   -- �����̩�٤�Ϳ������ --
    if (jadd==1) then
      n(x,z)=n1(z)*(1.0+jouran*jrate)+n0
    endif

!   -- �������٤�Ϳ������ --
    if (jadd==2) then
      z1=(dble(z)*dz-zpeak+zbot)/HO(z)+jouran*0.01
      n(x,z)=npeak*cos(chi)*exp(1.0-z1-exp(-z1)/cos(chi))+n0
    endif
      n00(x,z)=n(x,z)
  enddo
enddo

!-- ������������ --
do x=1,xmax
  do z=1,zmax
    if (exiswind==1) then
!if (z < zmax/4) then
!if (x > xmax*1/5) then
!      if (x < xmax*2/5) then
        U(x,z)=U0
!endif
!      else !if((z > 36) then
!        U(x,z)=-U0
!      else
!        U(x,z)=0.0
!      endif
!endif
    else
      U(x,z)=0.0
    endif
  enddo
enddo

!-- �ڥ���������ή������ --
do x=1,xmax
  do z=1,zmax
    if (exisped==1) then
!      if ((x-xmax/3.0)**2+(z-zmax/3.0)**2 < 9) then
        Eped(x,z)=E0(z)
!      else if((x-xmax*2.0/3.0)**2+(z-zmax/3.0)**2 < 9) then
!        Eped(x,z)=-E0
!      else
!        Eped(x,z)=0.0
!      endif
    else
      Eped(x,z)=0.0
    endif
  enddo
enddo

do x=1,xmax
  n(x,0) =n(x,1)**2/n(x,2)
  n(x,zmax+1) =n(x,zmax)**2/n(x,zmax-1)
enddo

do z=0,zmax+1
  n(0,z) =n(xmax,z)     !������
  n(xmax+1,z)=n(1,z)    ! 
enddo
do x=1,xmax
  vix(x,0)=mO*g/e/B
  viz(x,0)=mO*g/e/B
enddo
do z=0,zmax
  vix(0,z)=0
  viz(0,z)=0
enddo

!-- �ࡦE ������ --
do x=0,xmax+1
  do z=0,zmax+1
    divE(x,z)=0.0
  enddo
enddo

hassan=0

! writetest
open (22,file="n1.d")
do z=1,zmax
!write(*,*)z
  write(22,*) z,n1(z)
enddo
close(22)

!============================== ��������׻���ʬ =============================
open(14,file="/home/yumu/data/difn.d")
open(15,file="/home/yumu/data/bubble.d")
do t=1,tmax
bdet=0
sumdifn=0.0
!####--- v �η׻� ---####
  do x=1,xmax
    do z=1,zmax
      dnx=(n(x+1,z)-n(x-1,z))/2.0d0
      dnz=(n(x,z+1)-n(x,z-1))/2.0d0
      dphix=(phi(x+1,z)-phi(x-1,z))/2.0d0
      dphiz=(phi(x,z+1)-phi(x,z-1))/2.0d0
      vix(x,z)=mO*g/e/B &
&             -bi(z)/ki(z)**2*(dphix/dx+E0(z)) &
&             +alpha(z)*k*Tempi(z)/e/B/(n(x,z)+1.0e10)* &
& dnz/dz+1.0/(1.0+ki(z)**2)*U0
      viz(x,z)=-dphix/dx/B &
&             -bi(z)/ki(z)**2*(dphiz/dz+U0*B) &
&             -alpha(z)*k*Tempi(z)/e/B/(n(x,z)+1.0e10)*dnx/dx &
&             +ki(z)/(1.0+ki(z)**2)*U0
      vex(x,z)=0.0  !-U(x,z)
      vez(x,z)=-dphix/dx/B &
&             +alpha(z)*k*Tempi(z)/e/B/(n(x,z)+1.0e10)*dnx/dx 

    enddo
  enddo

! -- v �ζ����������� --
  do x=1,xmax 
    vix(x,0) =vix(x,1)
    vix(x,zmax+1) =vix(x,zmax)
    viz(x,0) =viz(x,1)
    viz(x,zmax+1) =viz(x,zmax)
    vex(x,0) =vex(x,1)
    vex(x,zmax+1) =vex(x,zmax)
    vez(x,0) =vez(x,1)
    vez(x,zmax+1) =vez(x,zmax)
  enddo
  do z=0,zmax+1
    vix(0,z) =vix(xmax,z)     ! �����������
    vix(xmax+1,z)=vix(1,z)    ! 
    viz(0,z) =viz(xmax,z)     ! �����������
    viz(xmax+1,z)=viz(1,z)    !
    vex(0,z) =vex(xmax,z)     ! �����������
    vex(xmax+1,z)=vex(1,z)    ! 
    vez(0,z) =vez(xmax,z)     ! �����������
    vez(xmax+1,z)=vez(1,z)    ! 
  enddo

  do x=1,xmax
    do z=0,zmax+1
      nvx(x,z)=n(x,z)*vix(x,z)
      nvz(x,z)=n(x,z)*viz(x,z)
    enddo
  enddo

!####--- n �η׻� ---####	
  do x=1,xmax
    do z=1,zmax

      if (vix(x,z) > 0) then
        dnx=n(x,z)-n(x-1,z)
      else
        dnx=n(x+1,z)-n(x,z)
      endif

      if (viz(x,z) > 0) then
        dnz=n(x,z)-n(x,z-1)
      else
        dnz=n(x,z+1)-n(x,z)
      endif

!      dphix=(phi(x+1,z)-phi(x-1,z))/2.0d0
!      if (dphix > 0) then
!        dnz=(n(x,z+1)-n(x,z-1))/2.0d0
!      else
!        dnz=n(x,z)-n(x,z-1)
!      endif
!
!      if (n(x,z) > 0) then
        dvix=vix(x,z)-vix(x-1,z)
        dviz=viz(x,z)-viz(x,z-1)
!      else
!        dvix=vix(x+1,z)-vix(x,z)
!        dviz=viz(x,z+1)-viz(x,z)
!      endif

      dnvx=nvx(x,z)-nvx(x-1,z)
      dnvz=nvz(x,z)-nvz(x,z-1)

!      dn=-dt*(dnvx/dx+dnvz/dz)

      dn=-dt*(vix(x,z)*dnx/dx+viz(x,z)*dnz/dz) !&
!&        -dt*n(x,z)*(dvix/dx+dviz/dz)
      n(x,z)=n(x,z)+dn
      difn=(n(x,z)-n00(x,z))**2/n00(x,z)**2
      sumdifn=sumdifn+difn    

      if (x==int(xmax/2)) then
        if (bdet==0) then
          if (n(x,z) > bpoint) then
            bdet=1
            write(15,*)t,zbot+dz*dble(z)
          endif
        endif
      endif
    enddo
  enddo

! -- n �ζ����������� --
  do x=1,xmax 
    n(x,0) =n(x,1)
    n(x,zmax+1) =n(x,zmax)
  enddo
  do z=0,zmax+1
    n(0,z) =n(xmax,z)     ! �����������
    n(xmax+1,z)=n(1,z)    ! 
  enddo


!####--- j �η׻� ---####
  do x=1,xmax
    do z=1,zmax
      jx(x,z)=n(x,z)*e*(vix(x,z)-vex(x,z))
      jz(x,z)=n(x,z)*e*(viz(x,z)-vez(x,z))

!write(*,*) 'jx(',x,z,')=',jx(x,z)
    enddo
  enddo


! -- j �ζ����������� --
  do x=1,xmax 
    jx(x,0) =jx(x,1)
    jx(x,zmax+1) =jx(x,zmax)
    jz(x,0) =jz(x,1)
    jz(x,zmax+1) =jz(x,zmax)
  enddo
  do z=0,zmax+1
    jx(0,z) =jx(xmax,z)     ! �����������
    jx(xmax+1,z)=jx(1,z)    ! 
    jz(0,z) =jz(xmax,z)     ! �����������
    jz(xmax+1,z)=jz(1,z)    ! 
  enddo

!####--- �ࡦJ �η׻� ---####
  do x=1,xmax
    do z=1,zmax
      djx=(jx(x+1,z)-jx(x-1,z))/2.0
      djz=(jz(x,z+1)-jz(x,z-1))/2.0
      divj(x,z)=0.0!djx/dx+djz/dz
    enddo
  enddo

!####--- pedersen ��ή��η׻� ---####
  do x=1,xmax
    do z=1,zmax
      dnx=(n(x+1,z)-n(x-1,z))/2.0
      pedersen(x,z)=dnx/dx*bi(z)*alpha(z)*Eped(x,z)
    enddo
  enddo

!####--- ��������η׻� ---####
  do x=1,xmax
    do z=1,zmax
      dnz=(n(x,z+1)-n(x,z-1))/2.0
      wind(x,z)=U(x,z)*B*(dnz/dz*bi(z)*alpha(z) &
&                          +n(x,z)*bidash(z)*alpha(z) &
&                          +n(x,z)*bi(z)*beta(z))
!write(*,*) 'wind(',x,z,')=',wind(x,z)
    enddo
  enddo

!####--- phi �η׻� ---####
  if (hassan<1) then

    do x=0,xmax+1
      do z=0,zmax+1
        phiold(x,z)=phi(x,z)
        phi(x,z)=0.0d0
      enddo
    enddo

    do i=1,imax
      ic=0 ! ��«��Ƚ���ͤ�����

      do x=1,xmax
        do z=1,zmax
          ni=n(x,z)
          dnx=(n(x+1,z)-n(x-1,z))/2.0d0
          dnz=(n(x,z+1)-n(x,z-1))/2.0d0
          dphix=(phi(x+1,z)-phi(x-1,z))/2.0d0
          dphiz=(phi(x,z+1)-phi(x,z-1))/2.0d0
          phisumx=phi(x+1,z)+phi(x-1,z)
          phisumz=phi(x,z+1)+phi(x,z-1)

           phinew(x,z)=dx**2*dz**2/2.0/(dx**2+dz**2)*( &
&                       phisumx/dx**2+phisumz/dz**2 &
&                       +1.0/(ni+1.0e10)*(dnx/dx*(dphix/dx-E0(z)) &
&                               +dnz/dz*(dphiz/dz-U0*B) &
&                               -ki(z)**2/bi(z)*dnx/dx*g/omega))

!          phinew(x,z)=((dz**2*dnx*dphix+dx**2*dnz*dphiz)/8.0/ni &
!&                   +(dz**2*phisumx+dx**2*phisumz)/2.0 &
!&                   -mO*omega*g/4.0/e/nu(z)/ni*dx*dz**2*dnx)/(dx**2+dz**2)

!          phinew(x,z)=dx**2*dz**2/2.0/(dx**2+dz**2)*( &
!&                       phisumx/dx**2+phisumz/dz**2 &
!&                       +ki(z)**2/ni/bi(z)*( &
!&                          dnx/dx*bi(z)/ki(z)**2*(dphix/dx-E0) &
!&                         +dnz/dz*bi(z)/ki(z)**2*dphiz/dz &
!&                         +ni*(alpha(z)*bidash(z)+beta(z)*bi(z))* &
!&                          (dphiz/dz-U0*B) &
!&                         -dnx/dx*mO*g/B/e))



!write(*,*) 'a=',dx**2*dz**2/2.0/(dx**2+dz**2)
!write(*,*) 'b=',ki(z)**2/ni/bi(z)
!write(*,*) 'c=',divj(x,z)/e
!         -- ��«��Ƚ�� --
          if (abs(phinew(x,z)-phi(x,z)) < abs(0.01*phi(x,z))) then
            ic=ic+1 
          endif

        enddo   ! z �Υ롼�׽�λ
      enddo   ! x �Υ롼�׽�λ

      do x=1,xmax
        do z=1,zmax
          phi(x,z)=phinew(x,z)
        enddo
      enddo

!     -- phi �ζ����������� --
      do x=1,xmax
        phi(x,0) =0.0 ! ��ü���Ű� = 0 �˸���
        phi(x,zmax+1) =0.0 !phi(x,zmax,t) !**2/phi(x,zmax-1,t)
      enddo
      do z=0,zmax+1  
        phi(0,z) =phi(zmax,z)
        phi(zmax+1,z)=phi(1,z)    
      enddo

!     -- ��«��Ƚ�� (���ʻ�) --
      if (ic==xmax*zmax) then
        goto 100
      endif

    enddo   ! i �Υ롼�׽�λ

    write(*,*) 'hassan'
    write(*,*) ic
!    hassan=hassan+1
    hassant=t

  else

    do x=0,xmax+1
      do z=0,zmax+1
        phi(x,z)=phiold(x,z)
      enddo
    enddo

  endif

100 continue

!if (t==???) then
! hassan=hassan+1
!endif

  write(*,*) 't=',t,'i=',i


!####--- ���� ---####

  if (mod(t,10)==0) then
    write(filename,'(a,i5.5,a)') '/home/yumu/data/rt', t, '.d'
    open(10,file=filename, status='unknown')
    write(filename2,'(a,i5.5,a)') '/home/yumu/data/phi', t, '.d' 
    open(11,file=filename2, status='unknown')
    write(filename3,'(a,i5.5,a)') '/home/yumu/data/vix', t, '.d' 
    open(12,file=filename3, status='unknown')
    write(filename4,'(a,i5.5,a)') '/home/yumu/data/viz', t, '.d' 
    open(13,file=filename4, status='unknown')

    do x=1,xmax
      do z=1,zmax
        write(10,*) (x*dx-xscale/2.0)/1.0e3, (z*dz+zbot)/1.0e3, n(x,z)
        write(11,*) (x*dx-xscale/2.0)/1.0e3, (z*dz+zbot)/1.0e3, phi(x,z)
!        write(12,*) (x*dx-xscale/2.0)/1.0e3, (z*dz+zbot)/1.0e3, vix(x,z)
!        write(13,*) (x*dx-xscale/2.0)/1.0e3, (z*dz+zbot)/1.0e3, viz(x,z)
      enddo 
      write(10,*) ! ����Ԥ�����
      write(11,*)
      write(12,*)
      write(13,*)
    enddo

    close(10)
    close(11)
    close(12)
  endif
write(14,*) t,sumdifn
enddo ! t �Υ롼�׽�λ
close(14)
close(15)
write(*,*) hassant
end
