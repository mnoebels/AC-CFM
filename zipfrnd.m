function out = zipfrnd(rho, sampleSize)
    a = rho;
    b = 2^(a-1);
    am1 = a - 1;
    bm1 = b - 1;
    u1 = rand( sampleSize );
    out = floor( u1.^(-1/am1) );
    clear('u1');
    u2 = rand( sampleSize );
    t = ( 1 + 1./out ).^(a-1);
    indxs = find( u2.*out.*(t-1)/bm1 > t/b );
    while ~isempty(indxs)
         indxsSize = size( indxs );
         u1 = rand( indxsSize );
         outNew = floor( u1.^(-1/am1) );
         clear('u1');
         u2 = rand( indxsSize );
         t = ( 1 + 1./outNew ).^(a-1);
         l = u2.*outNew.*(t-1)/bm1 <= t/b;
         out( indxs(l) ) = outNew(l);
         indxs = indxs(~l);
    end
end