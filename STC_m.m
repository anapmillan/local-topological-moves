function a = STC_m(N,q,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If you use this code, please cite:
%  A. P. Millan, R. Ghorbanchian, N. Defenu, F. Battiston and G. Bianconi
% "Local topological moves determine global diffusion properties of hyperbolic higher-order networks"
%  Physical Review E 104, 054302 (2021).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code that generates an STC network with m and q as given with a fixed size
% (number of nodes) equal to N.

% Output: a adjacency matrix
% Inputs: 
% Number of nodes: N
% Probability of closing a triangle: q
% Number of links per new node: m
% (c)  Ana P Millan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables
    imax = min(max(0.001*N,200),N); % imax is estimation for kmax.   
    
    %Neighbour matrix:
    v_neig = zeros(N,imax); 
    for i = 1:m+1
        ii = 1:m+1;
        ii(i) = [];
        v_neig(i,1:m) = ii;
        for j = i+1:m+1
            link_sat(i,j) = 1;
        end
    end

    % Degree array:
    ki = zeros(N,1);
    ki(1:m+1) = m;

    % Second neighbours:
    k2i = zeros(N,1); 
    v_neig2 = zeros(N,imax); 
 

%% Temporal loop: grow network
    for i = (m+2):N
        if mod(i,1e4)==0; fprintf('Node %d of %d\n',i,N); end

        % Pick first node to connect at random
        n1 = randi(i-1);

        % Pick process to select the other m-1 nodes
        aux = rand(m-1,1);
        wp = q > aux; %1 for process 1 (triangle), 0 for process 2 (square)
        wp2 = find(~wp);
        n_p2 = numel(wp2);

        %Make sure that we have enough 2nd neighbors, change to 1st otherwise
        while ( n_p2>k2i(n1) )
            aux = 1;
            wp(wp2(aux)) = 1;
            wp2 = find(~wp);
            n_p2 = numel(wp2);
        end

        %Select the 2nd node from 1st or 2nd neighbours according to process
        n2 = zeros(m-1,1);
        pos_nn1 = v_neig(n1,1:ki(n1)); %posible nodes for process 1: neighbours of n1
        pos_nn2 = v_neig2(n1,1:k2i(n1)); %posible nodes for process 2: 2nd neighbours of n1

        if sum(ismember(pos_nn1,pos_nn2))>0 %Sanity check
            error('Repeated!')
        end
        
        % Select m-1 nodes, each according to its process
        for j = 1:m-1
            if wp(j) %Process 1
                ival = randi(numel(pos_nn1)); %select index
                n2(j) = pos_nn1(ival); 
                pos_nn1(ival) = []; %Remove from list to avoid double pick

            else %Process 2
                ival = randi(numel(pos_nn2));
                n2(j) = pos_nn2(ival);
                pos_nn2(ival) = [];
            end
        end

        %% Update neighbour and second neigbour matrix and degrees

        % First neighbours
        %Node n1 
        v_neig(n1,ki(n1)+1) = i;
        ki(n1) = ki(n1)+1;

        %New node i
        v_neig(i,1:m) = sort([n1;n2]);
        ki(i) = m;

        %Nodes n2
        for j = 1:m-1
            v_neig(n2(j),ki(n2(j))+1) = i;
        end
        ki(n2) = ki(n2)+1;

        %2nd neighbors: neighbors of n1 and n2 (all)
        aux_v2 = v_neig(n1,1:ki(n1)-1);
        for j = 1:m-1
            aux_v2 = [aux_v2, v_neig(n2(j),1:ki(n2(j))-1)];
        end
        aux_v2 = unique(aux_v2);
        if aux_v2(1) == 0
            aux_v2(1) = [];
        end
        %Remove nodes that are already 1st neighbours
        iaux = ismember(aux_v2,[n1;n2]);
        aux_v2(iaux) = [];    
        %Add to new node
        v_neig2(i,1:numel(aux_v2)) = aux_v2;
        k2i(i) = numel(aux_v2);

        %Add new node to existing nodes 2nd neighbours
        for j = 1:numel(aux_v2) 
            v_neig2(aux_v2(j),k2i(aux_v2(j))+1) = i;
            k2i(aux_v2(j)) = k2i(aux_v2(j)) + 1;
        end

        %Neighbours of node i are now (at least) 2nn 
        for j = 1:ki(i) 
            aux = v_neig(i,1:ki(i));
            n_node = aux(j);
            aux(j) = []; 
            aux = aux(~ismember(aux,v_neig(n_node,1:ki(n_node)))); %Remove neighbours
            aux = unique([v_neig2(n_node,1:k2i(n_node)), aux]);  %Combine with old and Remove reps

            v_neig2(n_node,1:numel(aux)) = aux;
            k2i(n_node) = numel(aux);
        end

    end

    % Recover adjacency matrix
    a = zeros(N);
    for j = 1:N
        j2 = v_neig(j,1:ki(j));
        a(j,j2) =1;
    end
end
