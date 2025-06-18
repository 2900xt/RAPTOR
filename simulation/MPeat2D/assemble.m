function [ Kg ] = assemble( Kg, Kl, local_data, totalelementdof );
% assemble: To assemble the global stiffness matrix [Kg]  
%

     for ilocaldof=1:totalelementdof ;
         iglobaldof = local_data.dofs(1,ilocaldof);
        for jlocaldof=1:totalelementdof;
            jglobaldof = local_data.dofs(1,jlocaldof);
                Kg(iglobaldof,jglobaldof) =  Kg(iglobaldof,jglobaldof) + Kl(ilocaldof,jlocaldof);
        end
     end

