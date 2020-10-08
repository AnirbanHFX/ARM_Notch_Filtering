% Owner
% Anirban Majumder
% Git : https://github.com/AnirbanHFX
% Provided as is
%%
function out = PSO_fp(args)

    % Particle Swarm Optimization

    %%
    % Parameters of the solution space
    dim = args.dim;        % Dimension of the position vector
    order = args.order;    % Order of filter
    initrange.maxP = [];
    initrange.minP = [];
    initrange.maxV = [];
    initrange.minV = [];
    range = repmat(initrange, dim, 1);      % Structure of ranges of each position and velocity component
    for i=1:dim
        range(i).maxP = args.range(i).maxP;
        range(i).minP = args.range(i).minP;
        range(i).maxV = 0.3*(range(i).maxP - range(i).minP);
        range(i).minV = -1*range(i).maxV;
    end

    %%
    % Parameters of PSO
    iMax = args.iMax;            % Max iterations
    PopSize = args.PopSize;      % Number of particles
    w = args.w;                  % Inertia coefficient
    c1 = args.c1;                % Self best acceleration coefficient
    c2 = args.c2;                % Global best acceleration coefficient
    b = args.b;                  % Damping ratio

    % Variables for storing result of PSO
    global_costMin = inf;
    global_costMin_x = [];
    
    %%
    % Initialize particles
    p.x = zeros(order, dim);
    p.v = zeros(order, dim);
    p.cost = [];
    p.costMin = [];
    p.costMin_x = [];
    p_arr = repmat(p, PopSize, 1);

    for i=1:PopSize

        for k=1:order
            for j=1:dim
                 p_arr(i).x(k, j) = range(j).minP + (range(j).maxP - range(j).minP)*rand();
            end
        end

        p_arr(i).cost = costFn(p_arr(i).x, args);
        p_arr(i).costMin = p_arr(i).cost;
        p_arr(i).costMin_x = p_arr(i).x;

        if global_costMin > p_arr(i).costMin
            global_costMin = p_arr(i).costMin;
            global_costMin_x = p_arr(i).costMin_x;
        end

    end
    
    %%
    % Main Loop

    for idx=1:iMax

        for i=1:PopSize

            p_arr(i).v = w*p_arr(i).v + c1*rand(order,dim).*(p_arr(i).costMin_x - p_arr(i).x) + c2*rand(order,dim).*(global_costMin_x - p_arr(i).x);    % Update velocity
            
            for k=1:order       % Check for velocity saturation
                for j=1:dim
                    if p_arr(i).v(k, j) > range(j).maxV
                        p_arr(i).v(k, j) = range(j).maxV;
                    elseif p_arr(i).v(k, j) < range(j).minV
                        p_arr(i).v(k, j) = range(j).minV;
                    end
                end
            end

            p_arr(i).x = round(p_arr(i).x + p_arr(i).v);        % Update position

            for k=1:order       % Check for position bounds
                for j=1:dim
                    if p_arr(i).x(k, j) > range(j).maxP
                        p_arr(i).x(k, j) = range(j).maxP;
                    elseif p_arr(i).x(k, j) < range(j).minP
                        p_arr(i).x(k, j) = range(j).minP;
                    end
                end
            end

            p_arr(i).cost = costFn(p_arr(i).x, args);   % Compute cost

            if p_arr(i).cost < p_arr(i).costMin     % Update personal best

                p_arr(i).costMin = p_arr(i).cost;
                p_arr(i).costMin_x = p_arr(i).x;

                if p_arr(i).costMin < global_costMin    % Update global best

                    global_costMin = p_arr(i).costMin;
                    global_costMin_x = p_arr(i).costMin_x;

                end            

            end

        end

        w = w*b;

        disp(['Best cost = ' num2str(global_costMin)]);

    end

    %%
    % Return best X
    
    disp('Best X =')
    disp(global_costMin_x);
    
    X = zeros(args.order, 1);
    
    for i=1:args.order          % Reconstruct X from particle position
        X(i) = X(i) + global_costMin_x(i, 1);
        X(i) = X(i) + global_costMin_x(i, 2)/(2^args.fl);
    end
    ret = [];
    ret.X = X;
    ret.cost = global_costMin;
    out = ret;                  % Return best X
    
end



function y = costFn(x, args)    % Cost function: y = X^t*Q*X + P^t*X

    X = zeros(args.order, 1);
    
    for i=1:args.order          % Reconstruct column vector X from particle position
        X(i) = X(i) + x(i, 1);
        X(i) = X(i) + x(i, 2)/(2^args.fl);
    end
    
    y = double(transpose(X)*args.Q*X + transpose(args.P)*X);        % Return cost
    
end