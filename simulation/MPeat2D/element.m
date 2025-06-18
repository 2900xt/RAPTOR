function [ Kl,Cl,Ml ] = element( local_data )
% Element subroutine to construct local stiffness matrix from local data
%
  Kl = zeros(9,9);
  Cl = zeros(9,9);
  Ml = zeros(9,9);
  
  D  = zeros(3,3);
  A  = zeros(3,6);
  Df = zeros(2,2);
  Af = zeros(2,3);
%
%  Constitutive matrix - Hooke's Law and Darcy's Law
%
  Emod    = local_data.mater(1,1);
  nu      = local_data.mater(1,2);
  k       = local_data.mater(1,3);  
  sw      = local_data.mater(1,4);
  ss      = local_data.mater(1,5);
  poro    = local_data.mater(1,6);
%-------------------------------- Evalaute constitutive matrices
  D = sw*Emod/((1+nu)*(1-2*nu)) * [(1-nu)       nu      0;
                                  nu       (1-nu)    0;
                                  0          0     (1-2*nu)/2];
  Df= k* [ 1  0;
              0  1];
  m = [1;1;0];
  bint = (1/3)*[1 1 1];
%-------------------------------- Evaluate coefs and determinant for element area
    b = local_data.coords(2,1)*local_data.coords(3,2) - local_data.coords(3,1)*local_data.coords(2,2);
    b21 = local_data.coords(2,2) - local_data.coords(3,2);
    b22 = local_data.coords(3,2) - local_data.coords(1,2); 
    b23 = local_data.coords(1,2) - local_data.coords(2,2);
    b31 = local_data.coords(3,1) - local_data.coords(2,1);
    b32 = local_data.coords(1,1) - local_data.coords(3,1);
    b33 = local_data.coords(2,1) - local_data.coords(1,1);
    d2  = local_data.coords(1,1)*b21 + local_data.coords(1,2)*b31 + b;
    area = d2/2;
%-------------------------------- Evaluate derivatives of shape functions - mechanics
A = (1/d2) * [ b21   0   b22   0   b23   0;
                0   b31   0   b32   0   b33;
               b31  b21  b32  b22  b33  b23] ;
%-------------------------------- Evaluate derivatives of shape functions - flow
Af= (1/d2) * [ b21   b22   b23;
               b31   b32   b33];              
%-------------------------------- Assemble local stiffness matrix [ ux1 ux2 ux3 uy1 uy2 uy3 p1 p2 p3 ]
    thickness = 1;
    Kl = [ zeros(6,6)              zeros(6,3);
           zeros(3,6)              Af'*(Df*Af)*area*thickness  ];

%-------------------------------- Reorder stiffness matrix for dof

% [ux1 uy1 ux2 uy2 ux3 uy3 p1 p2 p3] ==> [ux1 uy1 p1 ux2 uy2 p2 ux3 uy3 p3]

% ----- Reorder columns  
            tempcol = Kl(:,3); Kl(:,3) = Kl(:,7); Kl(:,7) = tempcol;
            tempcol = Kl(:,6); Kl(:,6) = Kl(:,8); Kl(:,8) = tempcol;
            tempcol = Kl(:,4); Kl(:,4) = Kl(:,7); Kl(:,7) = tempcol;
            tempcol = Kl(:,5); Kl(:,5) = Kl(:,7); Kl(:,7) = tempcol;
% ----- Reorder rows
            temprow = Kl(3,:); Kl(3,:) = Kl(7,:); Kl(7,:) = temprow;
            temprow = Kl(6,:); Kl(6,:) = Kl(8,:); Kl(8,:) = temprow; 
            temprow = Kl(4,:); Kl(4,:) = Kl(7,:); Kl(7,:) = temprow;
            temprow = Kl(5,:); Kl(5,:) = Kl(7,:); Kl(7,:) = temprow;

%-------------------------------- Assemble local capacitance matrix [ ux1 ux2 ux3 uy1 uy2 uy3 p1 p2 p3 ]
    thickness = 1;
   B22 = (area*thickness*ss)*[1/6  1/12  1/12;
                            1/12  1/6  1/12;
                            1/12  1/12  1/6];
    Cl = [ A'*(D*A)*area*thickness               -A'*(m*bint)*area*thickness;
           -bint'*(m'*A)*area*thickness           B22   ];
           
%-------------------------------- Reorder capacitance matrix for dof

% [ux1 uy1 ux2 uy2 ux3 uy3 p1 p2 p3] ==> [ux1 uy1 p1 ux2 uy2 p2 ux3 uy3 p3]

% ----- Reorder columns  
            tempcol = Cl(:,3); Cl(:,3) = Cl(:,7); Cl(:,7) = tempcol;
            tempcol = Cl(:,6); Cl(:,6) = Cl(:,8); Cl(:,8) = tempcol;
            tempcol = Cl(:,4); Cl(:,4) = Cl(:,7); Cl(:,7) = tempcol;
            tempcol = Cl(:,5); Cl(:,5) = Cl(:,7); Cl(:,7) = tempcol;
% ----- Reorder rows
            temprow = Cl(3,:); Cl(3,:) = Cl(7,:); Cl(7,:) = temprow;
            temprow = Cl(6,:); Cl(6,:) = Cl(8,:); Cl(8,:) = temprow; 
            temprow = Cl(4,:); Cl(4,:) = Cl(7,:); Cl(7,:) = temprow;
            temprow = Cl(5,:); Cl(5,:) = Cl(7,:); Cl(7,:) = temprow;





         
         
