function k=wavek(T,h)

%k=wavek(T,h)
% computes radian wave number
% inputs are wave period and water depth
% T can be array & h is one depth
% h can be a vector & one wave period
% T and h can be equal arrays
% updated July 31, 2020 Jamie MacMahan

g=9.81; %gravity 

if length(T)==length(h)
    for  i=1:length(T) %period
        Ld=g.*T(i).^2./(2.*pi); K=2.*pi./Ld; sigma=2.*pi./T(i); 
        error=10;
        while error>0.000001
            error=sigma.^2-g.*K.*tanh(K.*h(i));
            K=K+0.000001;
        end
        k(i)=K;    
       
    end
elseif length(T)<length(h)
    for  i=1:length(h) 
        Ld=g.*T.^2./(2.*pi); K=2.*pi./Ld; sigma=2.*pi./T; 
        error=10;
        while error>0.00001
            error=sigma.^2-g.*K.*tanh(K.*h(i)); 
            K=K+0.00001;
        end
        k(i)=K;    
    end
elseif length(T)>length(h)   
    for i=1:length(T) 
        Ld=g.*T(i).^2./(2.*pi); K=2.*pi./Ld; sigma=2.*pi./T(i); 
        error=10;
        while error>0.000001 
            error=sigma.^2-g.*K.*tanh(K.*h); %inline function 
            K=K+0.000001;
        end
        k(i)=K;  
    end
end
