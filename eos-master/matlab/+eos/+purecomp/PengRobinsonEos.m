classdef PengRobinsonEos < eos.purecomp.CubicEosBase
    % PengRobinsonEos Peng-Robinson equation of state
    %
    %  This class provides methods to calculate thermodynamic properties
    %  based on Peng-Robinson equation of state.
    
    properties (Constant, Access = private)
        Sqrt2 = sqrt(2)
        Delta1 = 1 + sqrt(2)
        Delta2 = 1 - sqrt(2)
    end
    properties (SetAccess = private)
        AcentricFactor % Acentric factor
    end
    methods (Static)
        function coeffs = zFactorCubicEq(A,B)
            % Compute coefficients of Z-factor cubic equation
            %
            % coeffs = ZFACTORCUBICEQ(A,B)
            %
            % Parameters
            % ----------
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            %
            % Returns
            % -------
            % coeffs : Coefficients of the cubic equation of Z-factor
            arguments
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
            end
            coeffs = [1, B - 1, A - 2*B - 3*B^2, -A*B + B^2 + B^3];
        end
        function coeffs = dPdVPolyEq(T,a,b)
            % Compute coefficients of the polynomial of dPdV = 0.
            %
            % coeffs = DPDVPOLYEQ(T,a,b)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % a : Attraction parameter
            % b : Repulsion parameter
            %
            % Returns
            % -------
            % coeffs : Coefficients of the polynomial of dPdV = 0.
            arguments
                T (1,1) {mustBeNumeric}
                a (1,1) {mustBeNumeric}
                b (1,1) {mustBeNumeric}
            end
            R = eos.ThermodynamicConstants.Gas;
            coeffs = [R*T, 4*b*R*T - 2*a, 2*(b^2*R*T + a*b), ...
                2*b^2*(a - 2*b*R*T), b^3*(b*R*T - 2*a)];
        end
        function lnPhi = lnFugacityCoeffImpl(z,A,B)
            % Compute the natural log of fugacity coefficients.
            %
            % lnPhi = LNFUGACITYCOEFFIMPL(z,s)
            %
            % Parameters
            % ----------
            % z : Z-factor
            % A : Attraction parameter
            % B : Repulsion parameter
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
            arguments
                z (:,1) {mustBeNumeric}
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
            end
            Sqrt2 = eos.purecomp.PengRobinsonEos.Sqrt2;
            Delta1 = eos.purecomp.PengRobinsonEos.Delta1;
            Delta2 = eos.purecomp.PengRobinsonEos.Delta2;
            lnPhi = z - 1 - log(z - B) ...
                - A./(2*Sqrt2*B).*log((z + Delta1*B)./(z + Delta2*B));
        end
        function P = pressureImpl(T,V,a,b)
            % Compute pressure.
            %
            % P = PRESSUREIMPL(T,V,a,b)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % V : Volume [m3]
            % a : Attraction parameter
            % b : Repulsion parameter
            %
            % Returns
            % -------
            % P : Pressure [Pa]
            arguments
                T {mustBeNumeric}
                V {mustBeNumeric}
                a {mustBeNumeric}
                b {mustBeNumeric}
            end
            R = eos.ThermodynamicConstants.Gas;
            P = R*T./(V - b) - a./((V - b).*(V + b) + 2*b*V);
        end
    end
    methods
        function obj = PengRobinsonEos(Pc,Tc,omega,Mw)
            % Construct PR EoS
            %
            % obj = PENGROBINSONEOS(Pc,Tc,omega,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            %
            % Returns
            % -------
            % obj : PengRobinsonEos
            arguments
                Pc (1,1) {mustBeNumeric}
                Tc (1,1) {mustBeNumeric}
                omega (1,1) {mustBeNumeric}
                Mw (1,1) {mustBeNumeric}
            end
            obj@eos.purecomp.CubicEosBase(0.45724,0.07780,Pc,Tc,Mw);
            obj.AcentricFactor = omega;
        end
        function obj = setParams(obj,Pc,Tc,omega,Mw)
            % Set parameters.
            %
            % obj = obj.SETPARAMS(Pc,Tc,omega,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            %
            % Returns
            % -------
            % obj : PengRobinsonEos
            arguments
                obj {mustBeA(obj,'eos.purecomp.PengRobinsonEos')}
                Pc (1,1) {mustBeNumeric}
                Tc (1,1) {mustBeNumeric}
                omega (1,1) {mustBeNumeric}
                Mw (1,1) {mustBeNumeric}
            end
            obj = setParams@eos.purecomp.CubicEosBase(obj,Pc,Tc,Mw);
            obj.AcentricFactor = omega;
        end
        function alpha = temperatureCorrectionFactor(obj,Tr)
            % Compute temperature correction factor.
            %
            % alpha = obj.TEMPERATURECORRECTIONFACTOR(Tr)
            %
            % Parameters
            % ----------
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % alpha : Temperature correction factor
            arguments
                obj {mustBeA(obj,'eos.purecomp.PengRobinsonEos')}
                Tr (:,:) {mustBeNumeric}
            end
            omega = obj.AcentricFactor;
            m = 0.3796 + 1.485*omega - 0.1644*omega^2 + 0.01667*omega^3;
            alpha = (1 + m*(1 - sqrt(Tr))).^2;
        end
    end
end