integer function PBCSYMI(i,dimi) !returns the PBC cell coordinate
integer i, dimi, p
p = abs(i)+10
PBCSYMI = mod(i-1+p*dimi, dimi) + 1
end function

