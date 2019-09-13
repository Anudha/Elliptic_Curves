%This script generates Mordell Curves and enumerates (x,y) points (where x,y are integers) on the Mordell Curve
%Form of EC:  y^2 = x^3 +bx+c
%Subset, Form of Mordell Curves:  y^2 = x^3+n  , i.e. n<x or can also be written with bx term
%Anudha Mittal 

close all
clear all
clc
Z=50;
integers=1:Z; % use Z here and many other variables can be changed to Z, but wanted to keep it general incase not all integers are tested later
squares=integers.*integers;
cubes=integers.*squares;

n=NA(length(squares), length(cubes));
for i=1:length(squares)
  for j=1:length(cubes)
    n(i,j)=squares(i)-cubes(j);
  endfor
endfor
n;
n_vector=reshape(n,length(squares)*length(cubes),1);
range(n) ; % why are they all equal? It's the square of the max integer,Z-1
range(n_vector);
##figure
##hist(n_vector, range(n_vector))

[a,b] = histc(n_vector,unique(n_vector)); % a is a vector with count of instances for each unique value, b is an index
y = a(b); % y = number of times x repeats in vector; y has same dimensions as vector
y(n_vector==0)=0 ;%this line removes the count of n=0 
%NumberOfIntegerPointsOnCurve=4
NumberOfIntegerPointsOnCurve=max(y)
temp=find(y>NumberOfIntegerPointsOnCurve-1);
Desired_Mordell_n=n_vector(temp);
UniqueDesired_Mordell_n=unique(Desired_Mordell_n)
l=length(UniqueDesired_Mordell_n); %This can later be changed to a for loop and run through the length of UniqueDesired_Mordell_n
temp=find(n_vector==UniqueDesired_Mordell_n(l));

%temp=find(n_vector==225); %The value here comes from looking at above line, Remove the semicolon to get output of above line
for i=1:length(temp);
Mordell_n=n_vector(temp(i));
end
x_indices=mod(temp,Z);x_indices(x_indices==0)=Z;
y_indices=ceil(temp/Z);


Mordell_Curve_n=n(x_indices(1),y_indices(1))
check=squares(x_indices)-cubes(y_indices);
[x_ints]=integers(y_indices)
[y_ints]=integers(x_indices)
check_2=y_ints.^2-x_ints.^3;

%plot Mordell curve equations, ignore n=0 case
clear x,y;
clear n;
x=(-1-max(x_ints)):.1:(max(x_ints)+1);
n=Mordell_Curve_n;
y1=sqrt(x.^3+n);
y2=-y1;
figure
plot(x,y1,x,y2); hold on
plot(x_ints,y_ints,'o',x_ints,-y_ints,'o');

text1=strcat("n = ",num2str(n));
text2= strcat("(x,y) pairs on EC = ", num2str(NumberOfIntegerPointsOnCurve*2));
text3= strcat("For (x,y) pairs: 0 < x and x,y < ",num2str(Z));
title("Mordell Curve:  y^2 = x^3 + n")
text ((-max(x_ints)*0.8),(max(y1)*0.8),{text1;text2;text3})
xlabel('x')
ylabel('y')

%Next version: ignore n=0 case for large Z