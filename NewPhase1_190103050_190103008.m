clc
clear all

%input
func = input('Enter the function number (1,2,3,4,5,6) - ','s');
obj = input('Enter if you want to maximize(1) or minimize(2) the problem - ');
a = input('Enter the lower bound - ');
b = input('Enter the upper bound - ');
del = input('Enter the increment - ');
E = input('Enter the termination condition parameter - ');


%Assigning question
switch func 
    case '1'
        f= @f1;
    case '2'
        f= @f2;
    case '3'
        f= @f3;
    case '4'
        f= @f4;
    case '5'
        f= @f5;
    case '6'
        f= @f6;
    otherwise
        disp('Enter valid function!')
end

% Calling Bounding Phase Method
[afin,bfin,t1] = BP(f,obj,a,b,del);
disp("Upper and lower bounds are ("+afin+","+bfin+")")

% Calling Newton-Raphson Method
[min,t2] = NR(f,obj,afin,bfin,E);
disp("Minima is - "+min)

ques=input('Do you want to export the results? yes(1) or no(0): ');
if ques==1
    writetable(t1,'table.xlsx','Sheet',1);
    writetable(t2,'table.xlsx','Sheet',2);
end

% Bounding Phase Method Function
function [afin,bfin,t1] = BP(f,obj,a,b,del)

    while true
        
        %Step1
        x = rand(1,1)*(b-a) +a;
        
        %Step2
        if f(x-abs(del),obj)>=f(x,obj) && f(x,obj)>=f(x+abs(del),obj) 
            del = abs(del);
            break
        elseif f(x-abs(del),obj)<=f(x,obj) && f(x,obj)<=f(x+abs(del),obj) 
            del = -1*abs(del);
            break
        else 
            continue
        end
    end


    k = 0;
    kx = x;
    y=f(x,obj);
    
    %Step3
    while  true
        xnew = x + (2.^k)*del;
        y(end+1)=f(xnew,obj);

        
        if xnew<a
            xnew = a;
            kx(end+1) = xnew;
            
            afin = kx(end);
            bfin = kx(end-2);
            break

        elseif xnew>b
            xnew = b;
            kx(end+1) = xnew;
           
            afin = kx(end-2);
            bfin = kx(end);
            break
        
        %Step4
        elseif f(xnew,obj) < f(x,obj)
            k = k+1;
            x = xnew;
            kx(end+1) = xnew;
             
            continue
        else
            if del<0
                kx(end+1) = xnew;
              
                afin = kx(end);
                bfin = kx(end-2);
            elseif del>0
                kx(end+1) = xnew;
                
                afin = kx(end-2);
                bfin = kx(end);
            end
            
            break
        end
    end  
    plot(kx,y,'bo','DisplayName','Bounding Phase');
       
       xlabel('X')
       ylabel('Function evaluation')
         
      hold on

    Values_of_x= kx.';
    Function_Evaluation=y.';
    iter=1:size(Function_Evaluation);
    Iterations=iter.';
    t1=table(Iterations,Values_of_x,Function_Evaluation);

    
end

% Newton Raphson Method Function
function [min,t2] = NR(f,obj,afin,bfin,E)

    %Step1
    x = rand(1,1)*(bfin-afin) +afin;
    
    k=1;
    kx = x;
    y = [f(x,obj)];
    flow = f(x-E,obj);
    fhigh = f(x+E,obj);    
    F_dash = (fhigh -flow)/(2*E);
    %Step2
    while abs(F_dash)>E
       
        F_double_dash = (fhigh - 2*y(end) + flow)/(E*E);
    
        x = x - ((F_dash)/(F_double_dash));

        y(end+1)=f(x,obj);
    
        kx(end+1) = x;
       
        flow = f(x-E,obj);
        fhigh = f(x+E,obj);
    
        F_dash = (fhigh - flow)/(2*E);

        k = k+1;

  
    end

    min = kx(end);
    figure(1)
    plot(kx,y,'ro','DisplayName','Newton Raphson')
    hold on
    plot(min,y(end),'g*','DisplayName','Minima')
    hold off

    figure(2)
    l = 1:size(y.');
    plot(l,y,'b.','MarkerSize',10,'DisplayName','Function evaluations')
    xlabel('Iterations');
    ylabel('Function Evaluations')
    xlim([0 l(end)+1]);

    hold on
    plot(l(end),y(end),'r.','MarkerSize',10,'DisplayName','Minima')
    hold off
    
    legend ('show')
 
    Values_of_x= kx.';
    Function_Evaluation=y.';
    iter=1:size(Function_Evaluation);
    Iterations=iter.';
    t2=table(Iterations,Values_of_x,Function_Evaluation);
 
end

% Questions
function y = f1(x,obj)
   if obj == 2
       y = ((2*x-5).^4) - ((x.^2 - 1).^3);
   elseif obj == 1
       y = -1*(((2*x-5).^4) - ((x.^2 - 1).^3));

   end
   
end

function y = f2(x,obj)
    if obj == 2
        y = (8 + x.^3) - 2*x - 2*(exp(x));
    elseif obj == 1
        y = -1*((8 + x.^3) - 2*x - 2*(exp(x)));
    end
end

function y = f3(x,obj)
    if obj == 2
        y = 4*x*sin(x);
    elseif obj == 1
        y = -1*(4*x*sin(x));
    end
end

function y = f4(x,obj)
    if obj == 2
        y = 2*((x-3).^2) + (exp(0.5*x*x));
    elseif obj == 1
        y = -1*(2*((x-3).^2) + (exp(0.5*x*x)));
    end
end

function y = f5(x,obj)
    if obj == 2
        y = x*x - 10*(exp(0.1*x));
    elseif obj == 1
        y = -1*(x*x - 10*(exp(0.1*x)));
    end
end

function y = f6(x,obj)
    if obj == 2
        y = 20*sin(x) - 15*x*x;
    elseif obj == 1
        y = -1*(20*sin(x) - 15*x*x);
    end
end