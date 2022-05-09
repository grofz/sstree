  module point_mod
    implicit none
    private
    public point_t, POINT_DIM, WP, region_t, sphere_t, rectangle_t
    public point_mag

    integer, parameter :: WP = selected_real_kind(15, 307)
    integer, parameter :: POINT_DIM = 3

    ! Point class (public?)
    type point_t
      real(WP) :: x(POINT_DIM)
    contains
      procedure :: distance => point_distance
      procedure :: isin => point_isin
      procedure, private :: point_equals
      generic :: operator(==) => point_equals
    end type

    interface point_t
      module procedure new_point
    end interface

    interface point_mag
      module procedure point_mag2
      module procedure point_mag3
      module procedure point_mag4
    end interface



    ! region class and sub-classes
    type, abstract :: region_t
    contains
      procedure(region_intersectsPoint), deferred :: intersectsPoint
      procedure(region_intersectsRegion), deferred :: intersectsRegion
    end type

    abstract interface
      elemental function region_intersectsPoint(this, point) result(res)
        import
        logical :: res
        class(region_t), intent(in) :: this
        class(point_t), intent(in) :: point
      end function

      elemental function region_intersectsRegion(this, region) result(res)
        import
        logical :: res
        class(region_t), intent(in) :: this, region
      end function
    end interface



    type, extends(region_t) :: sphere_t
      type(point_t) :: center
      real(WP)      :: radius
    contains
      procedure :: intersectsPoint => sphere_intersectsPoint
      procedure :: intersectsRegion => sphere_intersectsRegion
    end type

    interface new_sphere
      module procedure :: new_sphere
    end interface



    type, extends(region_t) :: rectangle_t
      type(point_t) :: lcor, rcor
    contains
      procedure :: intersectsPoint => rectangle_intersectsPoint
      procedure :: intersectsRegion => rectangle_intersectsRegion
    end type

    interface rectangle_t
      module procedure new_rectangle
    end interface

  contains
!
! Class constructors
!
    pure type(point_t) function new_point(x)
      real(WP), intent(in) :: x(:)
      if (size(x) > POINT_DIM) &
          error stop 'new point - out of dimension bounds'
      new_point % x = 0.0
      new_point % x(1:size(x)) = x
    end function

    elemental type(sphere_t) function new_sphere(p,r)
      class(point_t), intent(in) :: p
      real(WP), intent(in) :: r
      new_sphere % center = p
      new_sphere % radius = r
    end function

    elemental type(rectangle_t) function new_rectangle(p0,p1)
      class(point_t), intent(in) :: p0, p1
      new_rectangle % lcor = p0
      new_rectangle % rcor = p1
    end function



!
! Point class methods
!
    elemental function point_distance(p0, p1) result(dis)
      real(WP) :: dis
      class(point_t), intent(in) :: p0, p1

      ! Euclidean distance between points P0 and P1
      ! The metric should follow these rules:
      ! (1) dis is non-negative
      ! (2) dis == 0 only for two same points
      ! (3) symmetry: dis(P0,P1) == dis(P1,P0) 
      ! (4) triangular inequality: dis(AB) <= dis(AC)+dis(CB),
      !     it equals only if C lies on the line AB 
      dis = sum((p0 % x - p1 % x)**2)
      dis = sqrt(dis)
    end function point_distance



    elemental function point_mag2(p0, p1) result(dis)
      real(WP) :: dis
      class(point_t), intent(in) :: p0, p1
!
! Length of line P0-P1
!
      dis = sum((p0 % x - p1 % x)**2)
      dis = sqrt(dis)
    end function point_mag2



    elemental function point_mag3(p0, p1, p2) result(area)
      real(WP) :: area
      class(point_t), intent(in) :: p0, p1, p2
!
! Area of triangle P0-P1-P2
!
      real(WP), dimension(POINT_DIM) :: e01, e02, xvec
      e01 = p1 % x - p0 % x
      e02 = p2 % x - p0 % x
      xvec = vec_product(e01, e02)
      area = sqrt(sum(xvec**2))*0.5_WP
    end function point_mag3



    elemental function point_mag4(p0, p1, p2, p3) result(volu)
      real(WP) :: volu
      class(point_t), intent(in) :: p0, p1, p2, p3
!
! Volume of tetrahedra P0/P1/P2/P3 (P3 must on right-hand side loop P0/P1/P2)
!
      real(WP), dimension(POINT_DIM) :: e01, e02, nvec
      real(WP) :: leng, area, pp, heig
      e01 = p1 % x - p0 % x
      e02 = p2 % x - p0 % x
      nvec = vec_product(e01, e02)
      leng = sqrt(sum(nvec**2))
      area = leng * 0.5_WP
      nvec = nvec / leng
      pp = -dot_product(nvec, p0 % x)
      heig = dot_product(nvec, p3 % x) + pp
      volu = area * heig / 3.0_WP
    end function point_mag4



    elemental logical function point_equals(p0, p1)
      class(point_t), intent(in) :: p0, p1
      point_equals = all(p0 % x == p1 % x)
    end function point_equals



    pure integer function point_isin(this, arr)
      class(point_t), intent(in) :: this, arr(:)

      integer :: i

      point_isin = -1
      do i = 1, size(arr)
        if (arr(i) == this) then
          point_isin = i
          exit
        endif
      enddo
    end function



    pure function vec_product(a, b) result(c)
      real(WP), intent(in) :: a(POINT_DIM), b(POINT_DIM)
      real(WP) :: c(POINT_DIM)
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
    end function vec_product




!
! intersectsPoint / intersectsRegion methods
!
    elemental function sphere_intersectsPoint(this, point) result(res)
      logical :: res
      class(sphere_t), intent(in) :: this
      class(point_t), intent(in) :: point

      ! point-sphere intersection?
      res = this % center % distance(point) < this % radius
    end function



    elemental function rectangle_intersectsPoint(this, point) result(res)
      logical :: res
      class(rectangle_t), intent(in) :: this
      class(point_t), intent(in) :: point

      ! point-rectangle intersection?
      associate (a => this % lcor % x, &
                 b => this % rcor % x, &
                 p => point % x)
        res = all(a <= p .and. p <= b)
      end associate
    end function



    elemental function sphere_intersectsRegion(this, region) result(res)
      logical :: res
      class(sphere_t), intent(in) :: this
      class(region_t), intent(in) :: region

      select type(region)
      class is (sphere_t)
        res = intersects_SS(this, region)
      class is (rectangle_t)
        res = intersects_RS(region, this)
      class default
        error stop 'sphere_intersectsRegion - uknown type'
      end select
    end function



    elemental function rectangle_intersectsRegion(this, region) result(res)
      logical :: res
      class(rectangle_t), intent(in) :: this
      class(region_t), intent(in) :: region

      select type(region)
      class is (sphere_t)
        res = intersects_RS(this, region)
      class is (rectangle_t)
        res = intersects_RR(this, region)
      class default
        error stop 'rectangle_intersectsRegion - uknown type'
      end select
    end function


!TODO test!??
    elemental function intersects_RS(r, s) result(res)
      logical :: res
      class(rectangle_t), intent(in) :: r
      class(sphere_t), intent(in) :: s

      ! rectangle-sphere intersection?
      type(point_t) :: wid ! rectangle size/2 vector
      type(point_t) :: cen ! center of rectangle (C.o.R.)
      type(point_t) :: r2s ! vector from C.o.R. to sphere center 
      real(WP) :: dis

      wid % x = abs(r % rcor % x - r % lcor % x) * 0.5_WP
      cen % x =    (r % lcor % x + r % rcor % x) * 0.5_WP
      r2s % x = abs(s % center % x - cen % x)

      ! "wid" and "r2s" have positive components 
      if (any(r2s % x > (wid % x + s % radius))) then
        res = .false.
      elseif (any(r2s % x <= wid % x)) then
        res = .true.
      else
        res = r2s % distance(wid) <= s % radius
      endif
    end function



    elemental function intersects_SS(s1, s2) result(res)
      logical :: res
      class(sphere_t), intent(in) :: s1, s2

      ! sphere-sphere intersection?
      res = s1 % center % distance(s2 % center) < (s1 % radius + s2 % radius)
    end function



    elemental function intersects_RR(r1, r2) result(res)
      logical :: res
      class(rectangle_t), intent(in) :: r1, r2

      ! rectangle-rectangle intersection?
      associate (a => r1 % lcor % x, &
                 b => r1 % rcor % x, &
                 c => r2 % lcor % x, &
                 d => r2 % rcor % x)
        ! Line ends sequences "A-B C-D" and "C-D A-B" are the only cases
        ! when segments AB and CD do not intersect. 
        ! Because "a<b" and "c<d" is always true, segments do not intersect
        ! iff "b<c" or "d<a". 
        res = all(.not.(b < c .or. d < a))
      end associate
    end function

  end module point_mod
