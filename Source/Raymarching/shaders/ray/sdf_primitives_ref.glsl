#version 330 core

//ship
const vec3 SHIP_CENTER = vec3(-5., 0., -2.);

//
////// This is just a list of different sdf and functions like soft_min and soft_max used to combine them together
//

mat3 Identity3 = mat3(1,0,0,
                            0,1,0,
                            0,0,1);

mat3 rot_x(float x) {
    return mat3(1,0,0,
                0,cos(x),-sin(x),
                0,sin(x), cos(x));
}

mat3 rot_y(float x) {
    return mat3(cos(x),0,sin(x),
                0,1,0,
                -sin(x),0, cos(x));
}

mat3 rot_z(float x) {
    return mat3(cos(x),-sin(x),0,
                sin(x),cos(x),0,
                0,0, 1);
}

mat3 rot_axis(vec3 u, float a) {
    return mat3(cos(a) + u.x*u.x*(1-cos(a)), u.x*u.y*(1-cos(a)) - u.z*sin(a), u.x*u.z*(1-cos(a)) + u.y*sin(a),
                u.y*u.x*(1-cos(a)) + u.z*sin(a), cos(a) + u.y*u.y*(1-cos(a)), u.y*u.z*(1-cos(a)) - u.x*sin(a),
                u.z*u.x*(1-cos(a)) - u.y*sin(a), u.z*u.y*(1-cos(a)) + u.x*sin(a), cos(a) + u.z*u.z*(1-cos(a)));
}

float sdfSphere(vec3 p, vec3 center, float r) {
  return length(p-center) - r;
}

float sdfBox(vec3 p, vec3 center, vec3 R, mat3 orientation) {
  vec3 q = abs(orientation*(p - center)) - R;

  float max_coord = max(max(q.x,q.y), q.z);

  return length(max(q, 0)) + min(max_coord, 0);
}

float sdfPlane(vec3 p, float height, mat3 orientation) {
  vec3 q = orientation*p;

  return q.y - height;
}

float sdfTorus( vec3 p, vec2 t, vec3 center, mat3 orientation)
{
  p = orientation*(p-center);
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdfCone( vec3 p, vec3 center, vec2 lh, mat3 orientation )
{
  p = orientation*(p-center);
  float hyp = length(lh);
  float l = lh.x;
  float h = lh.y;
  vec2 c = vec2(l/hyp, h/hyp);
  vec2 q = h*vec2(c.x/c.y,-1.0);
    
  vec2 w = vec2( length(p.xz), p.y );
  vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
  vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
  float k = sign( q.y );
  float d = min(dot( a, a ),dot(b, b));
  float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
  return sqrt(d)*sign(s);
}

float sdfCappedCylinder( vec3 p, float h, float r, vec3 center, mat3 orientation )
{
  p = orientation*(p-center);
  vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r,h);
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdfCone( vec3 p, vec3 center, vec2 lh, mat3 orientation )

float udfCone( vec3 p, vec3 c, float w, float h, mat3 orientation ) {
  vec3 q = orientation*(p-c);
  vec2 r = vec2(length(q.xz), q.y);
  float l = pow(h*h + w*w,0.5);
  vec2 tip = vec2(0,h);

  if (dot(r-tip,vec2(-w,h))>0) return length(r-tip);
  else return dot((r-tip), vec2(h/l, w/l));
}

float udfConeBound( vec3 p, vec3 c, float w, float h, mat3 orientation ) {
  vec3 q = orientation*(p-c);
  vec2 r = vec2(length(q.xz), q.y);
  vec2 s = vec2(length(q.xz), -q.y);

  if (r.y < 0) return length(max(s-vec2(w,0),0));
  else {
    float l = pow(h*h + w*w,0.5);
    vec2 tip = vec2(0,h);

    if (dot(r-tip,vec2(-w,h))>0) return length(r-tip);
    else return dot((r-tip), vec2(h/l, w/l));
  }
}

float soft_transition(float x, float y) {
    float k = 0.1;
    return max(k-abs(x-y), 0) * max(k-abs(x-y), 0) * max(k-abs(x-y), 0)/(6*k*k);

}

float soft_min(float x, float y) {
    return min(x,y) - soft_transition(x,y);
}

float soft_max(float x, float y) {
    return max(x,y) + soft_transition(x,y);
}

vec3 gradient(vec3 p) {
    vec3 dx = vec3(dX,0,0);
    vec3 dy = vec3(0,dX,0);
    vec3 dz = vec3(0,0,dX);
    vec3 grad;
    grad.x = (sdfSomething(p+dx)-sdfScene(p))/dX;
    grad.y = (sdfSomething(p+dy)-sdfScene(p))/dX;
    grad.z = (sdfSomething(p+dz)-sdfScene(p))/dX;
    return grad;
} 


//ship construction
float sdEllipsoid( vec3 p, vec3 r ) // approximated
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float sdCone( vec3 p, vec2 c, float h )
{
    vec2 q = h*vec2(c.x,-c.y)/c.y;
    vec2 w = vec2( length(p.xz), p.y );
    
	vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
    vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
    float k = sign( q.y );
    float d = min(dot( a, a ),dot(b, b));
    float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
	return sqrt(d)*sign(s);
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
  vec3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

float sdCappedCone(vec3 p, vec3 a, vec3 b, float ra, float rb)
{
    float rba  = rb-ra;
    float baba = dot(b-a,b-a);
    float papa = dot(p-a,p-a);
    float paba = dot(p-a,b-a)/baba;
    float x = sqrt( papa - paba*paba*baba );
    float cax = max(0.0,x-((paba<0.5)?ra:rb));
    float cay = abs(paba-0.5)-0.5;
    float k = rba*rba + baba;
    float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );
    float cbx = x-ra - f*rba;
    float cby = paba - f;
    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(cax*cax + cay*cay*baba,
                       cbx*cbx + cby*cby*baba) );
}


float sdRoundBox( vec3 p, vec3 b, float r )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
}


float opSubtract( float d1, float d2 ) {
	 return max(-d1,d2); 
}

vec3 opElongate(vec3 p, vec3 h )
{
    return p - clamp( p, -h, h );
}

float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); 

}


//give the sdf of the ship and the construction of it
float shipWindows(vec3 p) {
    return sdRoundBox(vec3(abs(p.x) - 0.56, p.y - 0.76, p.z-0.57), vec3(0.5, 0.5, 0.5), 0.05);
}

// returns 0 for body, 1 for window
float windowDeform(vec3 p){
	return 1.0 - smoothstep(0.0, 0.01, shipWindows(p));
}


float shipBody(vec3 p) {
	return opSmoothUnion(
		sdEllipsoid(p, vec3(0.3, 0.4, 0.3)),
		sdCapsule(vec3(abs(p.x), p.y, abs(p.z)), vec3(0.2, 0.0, 0.2), vec3(0.28, -0.1, 0.28), 0.015),
		0.02
	) + 0.01 * windowDeform(p);
}

float shipEngine(vec3 p){
	return opSubtract(
		sdCone(p + vec3(0.0, 0.34, 0.0), vec2(0.6, 1.0), 0.25),
		sdCone(p + vec3(0.0, 0.42, 0.0), vec2(0.40, 1.0), 0.05) - 0.07
	);
}


// ring and landing gear
float shipRing(vec3 p){
	return min(
		sdTorus(opElongate(p + vec3(0, 0.05, 0), vec3(0, 0.05, 0)), vec2(0.4, 0.013)),
		sdCappedCone(vec3(abs(p.x), p.y, abs(p.z)), vec3(0.1, -0.3, 0.1), vec3(0.17, -0.5, 0.17), 0.02, 0.01)
	);
}

float ship(vec3 pOrigin){
  vec3 p = pOrigin-SHIP_CENTER;
	return min(
		min(
			shipBody(p),
			shipEngine(p)
		),
		shipRing(p)
	);
}