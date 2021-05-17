classdef makeSpectrum
    % This class allows to directly evaluate the spectrum
    % I didn't go through a single call to an SDOF solver to make it faster. The constants are evaluated once and for all.
    properties
        signal
        time 
        dtSignal
        spectrum
        ariasFunction
        PGA
        envelopeShape
        damping=0.05;
    end
    
    properties(Constant)
        dtSpectrum=0.02;
        
    end
    properties(Dependent)
        period    
        A
        B       
    end
    
    methods
        function obj = makeSpectrum()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
                
        end
        function obj=set.time(obj, vecIn)
            obj.time=vecIn;
        end
        
%         function obj=get.spectrum(obj)
%             obj.spectrum;
%         end
        
        function obj=set.signal(obj, vecIn)
            obj.signal=vecIn;
        end
        
        function obj=set.damping(obj, temp)
            obj.damping=temp;
        end
        
        function obj=get.dtSignal(obj)
            obj=obj.time(2);
        end
        
        function arias=get.ariasFunction(obj)
            arias=(pi/(2*9.81))*obj.dtSignal*cumtrapz((obj.signal).^2);
        end
        
        function pga=get.PGA(obj)
            pga=max(abs(obj.signal));
        end
        
        function obj=get.period(obj)
            Tb=0.05; Tc=0.25; Td=1.2;
            n1=10;n=1;
            obj=[ (obj.dtSpectrum/n1:obj.dtSpectrum/n1:Tb) (Tb+obj.dtSpectrum/n1:obj.dtSpectrum/n1:Tc) (Tc+obj.dtSpectrum*n/2:obj.dtSpectrum*n/2:Td) (Td+obj.dtSpectrum*n:obj.dtSpectrum*n:10)]';
        end       
        
        function obj=get.envelopeShape(obj)
            % makes the accelerogram's envelope
            [pks,locs] = findpeaks(abs(obj.signal),obj.time,'MinPeakDistance',0.2);
            obj = spline(locs,pks,obj.time);
            obj=obj/max(obj);
        end
        
        function obj=get.A(obj)
            
            omg=2*pi./obj.period;
            omgd=omg.*(1-obj.damping^2)^0.5;
            S=sin(omgd*obj.dtSignal);
            C=cos(omgd*obj.dtSignal);
            O=exp(-obj.damping.*omg*obj.dtSignal);
            P=(1-obj.damping^2)^0.5;
            obj={O.*(obj.damping/P.*S+C) 1./omgd.*O.*S;-omg/P.*O.*S O.*(C-obj.damping/P.*S)};
        end
        
        function obj=get.B(obj)

            omg=2*pi./obj.period;
            omgd=omg.*(1-obj.damping^2)^0.5;
            S=sin(omgd*obj.dtSignal);
            C=cos(omgd*obj.dtSignal);
            O=exp(-obj.damping.*omg*obj.dtSignal);
            P=(1-obj.damping^2)^0.5;
            Q=(2*obj.damping^2-1)./(omg.^2*obj.dtSignal);
            R=2*obj.damping./(omg.^3*obj.dtSignal);            
            
            obj={O.*((Q+obj.damping./omg).*S./omgd+(R+1./omg.^2).*C)-R,...
                -O.*(Q.*S./omgd+R.*C)-1./omg.^2+R;...
                O.*((Q+obj.damping./omg).*(C-obj.damping./P.*S)-(R+1./omg.^2).*(omgd.*S+obj.damping.*omg.*C))+1./(omg.^2*obj.dtSignal),...
                -O.*(Q.*(C-obj.damping./P.*S)-R.*(omgd.*S+obj.damping.*omg.*C))-1./(omg.^2*obj.dtSignal)};
        end
        
        function obj=get.spectrum(obj)
            
            ag      = obj.signal;
            nPoint= length(ag);
            vec = zeros(1,length(obj.period))';
            omg=2*pi./obj.period;
            
            a11=obj.A{1,1};
            a12=obj.A{1,2};
            a21=obj.A{2,1};
            a22=obj.A{2,2};
            
            b11=obj.B{1,1};
            b12=obj.B{1,2};
            b21=obj.B{2,1};
            b22=obj.B{2,2};
            
            % I suggest you that implement a parfor here. I didn't just for sake of portability.
            for j=1:length(obj.period)
                
                d=zeros(1,nPoint);
                v=zeros(1,nPoint);                     
                
                f1=b11(j)*ag+ b12(j)*circshift(ag,1);
                f2=b21(j)*ag+ b22(j)*circshift(ag,1);
                % At this stage it could be implemented a matrix calculation for the two vectors.
                % However, to build the matrix takes up to 10 times the time required by the for loop.
                for i=1:nPoint-1
                    d(i+1)=a11(j)*d(i)+a12(j)*v(i)+ f1(i);
                    v(i+1)=a21(j)*d(i)+a22(j)*v(i)+ f2(i);
                    %    ga(i+1)=-(2*eta*omg(j)*v(i+1)+omg(j)^2*d(i+1));
                end
                
                ga=(2*obj.damping*omg(j).*v+omg(j)^2.*d);
                vec(j) = max(abs(ga));
                
            end
            vec(1)=max(abs(ag)) ;
            obj=vec;
        end
    end
end
