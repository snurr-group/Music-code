clear;

     % Given Isotherms calculates n/p and plots 

     %the file containing the variables
     load dat1
     load dat2
     
     %the number of data points in the file, dummy should be 8
     [length1,dummy]=size(dat1)
     [length2,dummy]=size(dat2)

     %the pressure and  sigma
     p1(1:length1) = dat1(1:length1,1); 
     p2(1:length2) = dat2(1:length2,1); 

     n1(1:length1) = dat1(1:length1,2); 
     n2(1:length2) = dat2(1:length2,2); 

     sig1(1:length1) = dat1(1:length1,3); 
     sig2(1:length2) = dat2(1:length2,3); 

     %n_p
     for i=2:length1;
	 n_p1(i)=n1(i)/p1(i);
     end
     n_p1(1)=n_p1(2);
     for i=2:length2;
	 n_p2(i)=n2(i)/p2(i);
     end
     n_p2(1)=n_p2(2);	 
     
     % the strings for graph points/lines
     s_1    = 'r-' ;
     s_2    = 'bd' ;
     plot(p1,sig1,s_1, p2, sig2, s_2);

     leg_1='pxylene';
     leg_2='mxylene';
     
     legend(leg_1, leg_2,0)
       
xlab='Pressure kpa';
xlabel(xlab);

ylab='\sigma - molecs/uc/kpa';
ylabel(ylab);

title('\sigma for pxylene-mxylene');

