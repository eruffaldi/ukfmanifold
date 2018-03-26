# TODO
# Support Projection
# Support SPD
# Support transport
import json
import sys

def makeop2(fx,sourcerange,sourcerange2,destrange):
	x1 = "x1({start}:{end})".format(start=sourcerange[0],end=sourcerange[1])
	x2 = "x2({start}:{end})".format(start=sourcerange2[0],end=sourcerange2[1])

	if type(fx) is str:
		content= fx.format(x1=x1,x2=x2)
	else:
		content = fx(x1,x2)
	return "\t\t\tz({start}:{end}) = {content};\n".format(start=destrange[0],end=destrange[1],content=content)

def makeop1(fx,sourcerange,destrange):
	x = "x1({start}:{end})".format(start=sourcerange[0],end=sourcerange[1])
	if type(fx) is str:
		content= fx.format(x=x)
	else:
		content = fx(x)
	return "\t\t\tz({start}:{end}) = {content};\n".format(start=destrange[0],end=destrange[1],content=content)

rotextra="""
function omega = so3log(R)

    theta = acos((trace(R)-1)/2); %acos(max(-1,min((trace(R)-1)/2,1)));
    if isreal(theta) == 0
        R = R/abs(det(R));
        theta =  acos((trace(R)-1)/2);
    end

    if abs(theta) < 1e-10
        B = 0.5;
        SO = (1/(2))*(R-R');  % =skew(omega)
        iV = eye(3); % use directly skew of omega
    else
        A = sin(theta)/theta;
        B = (1-cos(theta))/(theta*theta);
        SO = (1/(2*A))*(R-R');  % =skew(omega)
        %??
        % syms x real
        % A = sin(x)/x
        % B= (1-cos(x))/(x*x)
        % Q=1/(x*x)*(1 - A/2/B)
        % taylor(Q,x,0)
        %       x^4/30240 + x^2/720 + 1/12
        Q= 1/(theta^2)*(1 - A/2/B);
        iV = eye(3) - 1/2*SO + Q*SO^2; % use directly skew of omega
    end

    omega = [SO(3,2) SO(1,3) SO(2,1)];

end

% Emanuele Ruffaldi 2017 @ SSSA
function R = so3exp(omega)

theta = norm(omega);
if theta < 1e-12
    % Original
    %     A = 1;
    %     B = 0.5;
    %     C = 1/6;
    %     S = zeros(3);
    %     R = eye(3) + A*S + B*S^2;
    %     V = eye(3) + B*S + C*S^2;
    
        N = 10;
        R = eye(3);
        xM = eye(3);
        cmPhi = skew(omega);
        for n = 1:N
            xM = xM * (cmPhi / n);
            R = R + xM;
        end
        
        % Project the resulting rotation matrix back onto SO(3)
        R = R / sqrtm( R'*R ) ;
    
else
    %Original
    A = sin(theta)/theta;
    B = (1-cos(theta))/(theta^2);
    C = (theta-sin(theta))/(theta^3);
    S = skew(omega);
    R = eye(3) + A*S + B*S^2;
    %V = eye(3) + B*S + C*S^2;
    
    %Barfoot
    if 0==1    
        axis = omega/theta;
        cp = cos(theta);
        sp = sin(theta);
        sa = skew(axis);
        
        R = cp*eye(3) + (1-cp)*axis*(axis') + sp*sa;
    end
    assert(abs(det(R)-1) < 1e-5,'unitary');
end
end
			"""
se3matextra="""

function y = se3mat_inv(x)

	M = reshape(x,4,4);
	R = M(1:3,1:3);
	y = eye(4);
	y(1:3,1:3) = R';
	y(1:3,4) = -y(1:3,1:3)*M(1:3,4);
	y = reshape(y,1,[]);
end


function y = se3mat_log(x)
	
	M = reshape(x,4,4);
	t = M(1:3,4);
	R = M(1:3,1:3);

    theta = acos((trace(R)-1)/2); %acos(max(-1,min((trace(R)-1)/2,1)));
    if isreal(theta) == 0
        R = R/abs(det(R));
        theta =  acos((trace(R)-1)/2);
    end

    if abs(theta) < 1e-10
        B = 0.5;
        SO = (1/(2))*(R-R');  % =skew(omega)
        iV = eye(3); % use directly skew of omega
    else
        A = sin(theta)/theta;
        B = (1-cos(theta))/(theta*theta);
        SO = (1/(2*A))*(R-R');  % =skew(omega)
        %??
        % syms x real
        % A = sin(x)/x
        % B= (1-cos(x))/(x*x)
        % Q=1/(x*x)*(1 - A/2/B)
        % taylor(Q,x,0)
        %       x^4/30240 + x^2/720 + 1/12
        Q= 1/(theta^2)*(1 - A/2/B);
        iV = eye(3) - 1/2*SO + Q*SO^2; % use directly skew of omega
    end
    omega = [SO(3,2) SO(1,3) SO(2,1)];

	%y = [m1.log(getrot(x)), getpos(x)]; % not exact
	u = iV*t;
	y = [omega(:);u(:)]';
end

function y = se3mat_exp(x)

	omega = x(1:3);
	u = x(4:6); 

	theta = norm(omega);
	if theta < 1e-12
	    % Original
	    %     A = 1;
	    %     B = 0.5;
	    %     C = 1/6;
	    %     S = zeros(3);
	    %     R = eye(3) + A*S + B*S^2;
	    %     V = eye(3) + B*S + C*S^2;
	    
	        N = 10;
	        R = eye(3);
	        xM = eye(3);
	        cmPhi = skew(omega);
	        V = eye(3);
	        pxn = eye(3);
	        for n = 1:N
	            xM = xM * (cmPhi / n);
	            pxn = pxn*cmPhi/(n + 1);
	            R = R + xM;
	            V = V + pxn;
	        end
	        
	        % Project the resulting rotation matrix back onto SO(3)
	        R = R / sqrtm( R'*R ) ;
	    
	else
	    %Original
	    if 1==1
	        A = sin(theta)/theta;
	        B = (1-cos(theta))/(theta^2);
	        C = (theta-sin(theta))/(theta^3);
	        S = skew(omega);
	        R = eye(3) + A*S + B*S^2;
	        V = eye(3) + B*S + C*S^2;
	    else
	    %Barfoot
	        
	        axis = omega/theta;
	        cp = cos(theta);
	        sp = sin(theta);
	        cph = (1 - cos(theta))/theta;
	        sph = sin(theta)/theta;
	        sa = skew(axis);
	        
	        R = cp*eye(3) + (1-cp)*axis*(axis') + sp*sa;
	        V = sph * eye(3) + (1 - sph) * axis * (axis') + cph * sa;
	    end
	    
	end

	y = eye(4);
	y(1:3,1:3) = R;
	y(1:3,4) = V*u(:);
	y = reshape(y,1,[]);

	end

"""

se3quatextra=""

def solvemodel(m):
	if m["type"][0] == "Quat":
		return dict(product="qmult({x1},{x2})",delta="qmdelta({x1},{x2})",step="qmstep({x1},{x2})",inv="qconj({x})",pack="{x}",unpack="{x}")
	elif m["type"][0] == "Rot":
		return dict(extra=rotextra,product="reshape(reshape({x1},3,3)*reshape({x2},3,3),1,[])",delta="so3log(reshape({x1},3,3)*reshape({x2},3,3)')",step="reshape(so3exp({x2})*reshape({x1},3,3),1,[])",inv="reshape(reshap({x},3,3)',1,[])",exp="reshape(so3exp({x}),1,[])",log="so3log(reshape({x},3,3))",pack="reshape({x},1,[])",unpack="reshape({x},3,3)")
	elif m["type"][0] == "SE3Mat" or m["type"][0] == "SE3MatGlobal":
		return dict(extra=se3matextra,product="reshape(reshape({x1},4,4)*reshape({x2},4,4),1,[])",delta="se3mat_log(reshape({x1},4,4)*reshape(se3mat_inv({x2}),4,4))",step="reshape(reshape(se3mat_exp({x2}),4,4)*reshape({x1},4,4),1,[])",inv="se3mat_inv({x})",log="se3mat_log({x})",exp="se3mat_exp({x})",pack="reshape({x},1,[])",unpack="reshape({x},4,4)")
	#elif m["type"][0] == "SE3Quat":
	#	return dict(extra=se3quatextra,product="se3quat_prod({x1},{x2})",delta="se3log(se3mat_munflat({x1})*se3mat_munflat(se3inv({x2})))",step="se3mat_mflat(se3mat_munflat(m.exp({x2}))*se3mat_munflat({x1}))",inv="se3mat_inv({x})",log="se3mat_log({x})",exp="se3mat_exp({x})",pack="se3mat_mflat({x1})",unpack="se3mat_munflat({x1})")
	elif m["type"][0] == "Rn":
		n = m["type"][1]
		return dict(product="{x1}+{x2}",delta="{x1}-{x2}",step="{x1}+{x2}",inv="-{x}",log="{x}",exp="{x}",pack=lambda x: x,unpack=lambda x:x)
	else:
		raise Exception("!!! Unsupported %s %s" % (m["type"],m))
		#return dict(product="anyprod",delta="anydelta",step="anystep",inv="anyinv",log="anylog",exp="anyexp",pack=lambda x: x,unpack=lambda x:x)
def main():
	if len(sys.argv) != 2:
		print "expected input JSON";
		return
	jq = json.load(open(sys.argv[1],"rb"))
	output = jq.get("output","-")
	group = jq["group"]
	alg = jq["alg"]
	count = jq["count"]
	models = jq["models"]
	ginc = jq["groupinc"]
	ainc = jq["alginc"]
	print jq["type"],group,alg,models,jq["groupinc"],jq["alginc"]

	# check if 
	maybeexp="""m.log = @mlog;
				m.exp = @mexp;
			""" if jq.get("withlog",False) else ""

	solvedmodels = [solvemodel(m) for m in jq["models"]]
	output = """
		function m = customManifold()
			m = [];
			m.type = 'special';
			m.inv = @minvert;
			m.prod = @mproduct;
			{maybeexp}
			m.delta = @mdelta;
			m.step = @mstep;
			m.group = {group};
			m.alg = {alg};
			m.count = {count};
			m.pack = @mpack;
			m.unpack = @munpack;
            m.s = int_manisetup([],[],m);
		end
		""".format(group=group,maybeexp=maybeexp,alg=alg,count=jq["count"])

	# output inv: group -> group
	output += """
		function z = minvert(x1)
			z = zeros({group},1);
""".format(group=group)
	for i in range(0,len(jq["models"])):
		output += makeop1(solvedmodels[i]["inv"],(ginc[i]+1,ginc[i+1]),(ginc[i]+1,ginc[i+1]))
	output += """		end
"""
	# output prod: group,group -> group
	output += """
		function z = mproduct(x1,x2)
			z = zeros({group},1);
""".format(group=group)
	for i in range(0,len(jq["models"])):
		output += makeop2(solvedmodels[i]["product"],(ginc[i]+1,ginc[i+1]),(ginc[i]+1,ginc[i+1]),(ginc[i]+1,ginc[i+1]))
	output += """
		end
"""
	# output log: group -> alg
	if maybeexp != "":
		output += """
			function z = mlog(x1)
				z = zeros({alg},1);
	""".format(alg=alg)
		for i in range(0,len(jq["models"])):
			output += makeop1(solvedmodels[i]["log"],(ginc[i]+1,ginc[i+1]),(ainc[i]+1,ainc[i+1]))
		output += """
			end
	"""
	# output exp: alg -> group
	if maybeexp != "":
		output += """
			function z = mexp(x1)
				z = zeros({group},1);
	""".format(alg=alg,group=group)
		for i in range(0,len(jq["models"])):
			output += makeop1(solvedmodels[i]["exp"],(ainc[i]+1,ainc[i+1]),(ginc[i]+1,ginc[i+1]))
		output += """		end
	"""
	# output delta: group,group -> alg
	output += """
		function z = mdelta(x1,x2)
			z = zeros({alg},1);
""".format(alg=alg)
	for i in range(0,len(jq["models"])):
		output += makeop2(solvedmodels[i]["delta"],(ginc[i]+1,ginc[i+1]),(ginc[i]+1,ginc[i+1]),(ainc[i]+1,ainc[i+1]))
	output += """		end
"""
	# output step: group,alg -> group
	output += """
		function z = mstep(x1,x2)
			z = zeros({group},1);
""".format(group=group)
	for i in range(0,len(jq["models"])):
		output += makeop2(solvedmodels[i]["step"],(ginc[i]+1,ginc[i+1]),(ainc[i]+1,ainc[i+1]),(ginc[i]+1,ginc[i+1]))
	output += """		end
"""
	# output pack: cell array -> group
	output += """
		function z = mpack(x1)
			z = zeros({group},size(x1,2));
			for j=1:size(x1,2)
""".format(group=group)
	for i in range(0,len(jq["models"])):
		x = "x1{pre}{i},j{post}".format(pre="{",post="}",i=i+1)
		fx = solvedmodels[i]["pack"]
		if type(fx) is str:
			content = fx.format(x=x)
		else:
			content = fx(x)
		output += "\t\t\t\tz({start}:{end},j) = {content};\n".format(astart=ainc[i]+1,aend=ainc[i+1],start=ginc[i]+1,end=ginc[i+1],content=content)

	output += """			end
		end
"""
	# output unpack: group -> cell
	output += """
		function z = munpack(x1)
			z = cell({count},size(x1,2));
			for j=1:size(x1,2)
""".format(group=group,count=count)
	for i in range(0,len(jq["models"])):
		x = "x1({start}:{end},j)".format(start=ginc[i]+1,end=ginc[i+1])
		fx = solvedmodels[i]["unpack"]
		if type(fx) is str:
			content = fx.format(x=x)
		else:
			content = fx(x)
		output += "\t\t\t\tz{pre}{i},j{post} = {content};\n".format(pre="{",post="}",i=i+1,astart=ainc[i]+1,aend=ainc[i+1],start=ginc[i]+1,end=ginc[i+1],content=content)
	output += """			end
		end
"""

	ee = set()
	for m in solvedmodels:
		q= m.get("extra","")
		ee.add(q)
	for q in ee:
		output += q
		output += "\n"
	if output == "-":
		print output
	else:
		open(jq["output"],"w").write(output)

if __name__ == '__main__':
 	main()