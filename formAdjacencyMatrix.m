% function A = formAdjacencyMatrix(gcoord,node)


function A = formAdjacencyMatrix(gcoord,node)

[nnode,~] = size(gcoord);
[nel,nnel] = size(node);

num_entries = 2*nnel*nel;

i_v = zeros(num_entries,1); j_v = zeros(num_entries,1);
s_v = ones(num_entries,1);

ind = 1:nnel;
ind_n = circshift(ind,[0 -1]);

start_ind = 1; end_ind = nel;
for i = 1:nnel
    if(i>1)
        start_ind = start_ind + nel;
        end_ind = end_ind + nel;
    end
    
    i_v(start_ind:end_ind) = node(:,ind(i));
    j_v(start_ind:end_ind) = node(:,ind_n(i));
    start_ind = start_ind + nel;
    end_ind = end_ind + nel;
    i_v(start_ind:end_ind) = node(:,ind_n(i));
    j_v(start_ind:end_ind) = node(:,ind(i));
    
    
end


A = sparse(i_v,j_v,s_v,nnode,nnode);

A(A>1.0)=1.0;
