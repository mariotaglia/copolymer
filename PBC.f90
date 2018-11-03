integer function PBCSYMI(i,dimi) !returns the PBC cell coordinate
integer i, dimi
PBCSYMI = mod(i-1+5*dimi, dimi) + 1
end function

