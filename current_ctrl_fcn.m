% By: ZTG 2019-7-22

function [I,Int] = current_ctrl_fcn(V,V_ref,Kp,C_rate,capacity,dt,Int)
    I_chrgmax = -capacity*C_rate/1000; %A

    %% Original P-controller only
% 
%     % Compute control law current, but regulate if above max charge rate
%     I_ctrl = Kp*;
% 
%     if I_ctrl < I_chrgmax % charging current is negative, so higher charging rate is a more negative value
%         scale_Kp = I_chrgmax/I_ctrl; % bound the max current; here use max instead of min because charge current is (-)
%         % recompute Kp if         
%         I = scale_Kp*Kp*(V-V_ref);
%     else % case where charging current from controller doesn't exceed max rate
%         I = I_ctrl;
%     end
    
    
    %% PI implementation    
    Ki = 0.16;
    Kt = 0.1;
    
    %%% Input saturation reset form
    % Proportional & Derivative terms
    Prop = Kp*(V-V_ref);
    Deriv = 0;

    %Compute Current
    % Positive error between Vmax and V should correspond to a negative
    % (i.e. charging) current
    % I_ctrl = controller input to actuator; I(k) = actuator output
    I_ctrl = Prop + Int + Deriv;
    I = max(I_ctrl,I_chrgmax); % bound the max current; here use max instead of min because charge current is (-)
    
    % Update Integral Term:
    % controller-actuator error: 0 when no saturation
    err_sat = I - I_ctrl;  
    Int = Int + (Ki*(V-V_ref) + Kt*err_sat)*dt;
end