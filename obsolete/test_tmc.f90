
program main
implicit none

double precision, external :: TMC
integer :: nn, mm, n, m, n1, m1, n2, m2
double precision :: x1, x2

nn = 5
mm = 4
n = 3
m = -2

n1 = 6
m1 = -3
n2 = 1
m2 = 5

x1 = 0.1d0
x2 = 0.9d0

print '("<",4(i4,","),f12.6"|",4(i4,","),f12.6">")', nn, mm, n, m, x1, n1, m1, n2, m2, x2
print *, 2*(nn+n)+abs(mm)+abs(m),2*(n1+n2)+abs(m1)+abs(m2), mm+m, m1+m2
print *, TMC(nn, mm, n, m, n1, m1, n2, m2, sqrt(x2/x1))

end program
