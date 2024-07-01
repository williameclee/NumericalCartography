%% LIGHTSOURCE
% defines a light source
%
% Syntax
%   light = LightSource(type, position)
%   light = LightSource(type, position, color)
%
% Description
%   light = LightSource(type, position)
%       generates a light source at the position and with the specified 
%       type.
%   light = LightSource(type, position, color)
%       generates a light source at the position with the specified colour.
%
% Input Arguments
%   type - Type of light source
%       The options are 'parallel' and 'point'.
%   position - Position of the light source
%       - 1-by-2 vector [azimuth, altitude] in degrees for a parallel 
%           light.
%       - 1-by-3 vector [Px, Py, Pz] of the position.
%   color - Colour of the light source
%       A 1-by-3 vector [R, G, B] of the colour. The values are in the 
%       range of [0, 1].
%       The default value is [1, 1, 1].
%
% Output Arguments
%   light - Light source object
%
% Authored by
%   En-Chi Lee <williameclee@gmail.com>, 2024-07-01

classdef LightSource

    properties
        % Type of light source
        Type string {mustBeMember(Type, {'parallel', 'point'})};

        % Position of the light source
        Position(1, 3) double {mustBeNumeric};
        Direction(1, 3) double {mustBeNumeric};
        HorizontalDirection(1, 2) double {mustBeNumeric};
        
        % Azimuth and altitude of the light source
        Azimuth double {mustBeNumeric};
        Altitude double {mustBeNumeric};

        % Colour of the light source
        Color(1, 3) double {mustBeNumeric};
    end

    methods

        function obj = LightSource(type, varargin)
            %% Initialisation
            % Parse inputs
            p = inputParser;
            addRequired(p, 'Type', @(x) isstring(x) || ischar(x));
            addOptional(p, 'Position', [-1, 1, 0.5], ...
                @(x) isvector(x) && (length(x) == 2 || length(x) == 3));
            addOptional(p, 'Color', [1, 1, 1], ...
                @(x) isvector(x) && length(x) == 3);
            parse(p, type, varargin{:});
            type = p.Results.Type;
            position = p.Results.Position;
            color = p.Results.Color;

            %% Generating the LightSource
            % Assign the property
            obj.Type = string(lower(type));

            % Assign the light source
            if length(position) == 2 % Azimuth and altitude
                obj.Azimuth = position(1);
                obj.Altitude = position(2);

                % Update the light to be the vector form
                position = azald2vec(obj.Azimuth, obj.Altitude);
            else % Vector form
                [obj.Azimuth, obj.Altitude] = vec2azald(position);
            end

            switch obj.Type
                case 'point'
                    obj.Position = position;
                    obj.Direction = position / vecnorm(position);
                case 'parallel'
                    obj.Position = position / vecnorm(position);
                    obj.Direction = obj.Position;
            end

            % Assign the horizontal direction
            obj.HorizontalDirection = [obj.Direction(1), obj.Direction(2)];
            obj.HorizontalDirection = ...
                obj.HorizontalDirection / vecnorm(obj.HorizontalDirection);

            % Assign the colour
            obj.Color = color;

        end

    end

end
