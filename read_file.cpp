#include <iostream>


#include <osmium/handler.hpp>
#include <osmium/osm/node.hpp>
#include <osmium/osm/way.hpp>
#include <osmium/io/any_input.hpp>
#include <osmium/visitor.hpp>
#include <osmium/index/map/sparse_mem_array.hpp>
#include <osmium/handler/node_locations_for_ways.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>


#include  <map>

#include  "cgal_stuff.h"


// CAiro stuff
#include  <cairo/cairo.h>
#include  <cairo/cairo-pdf.h>
#include  <cairo/cairo-svg.h>



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#define  LONG_TO_METERS   111321.0


// Edge weight.
typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
typedef boost::property<boost::vertex_distance_t, double> VertexDistProperty;

//struct edge_properties
//{
//    double      weight ;
//};



// Graph.
typedef boost::adjacency_list< boost::listS,
                               boost::vecS,
                               boost::undirectedS,
                               VertexDistProperty,
                               //boost::no_property,
                               //EdgeProperties
                               EdgeWeightProperty
                               > Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::edge_descriptor  Edge;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;



using  vertices_id_map_t = std::map <long, int>;



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


class  Mapper
{
private:
    Box2d  bb;
    double   xscale, yscale;
    double   frame;
    double  width, height;
    
public:
    Mapper( Box2d  &_bb ) {
        bb = _bb;
        yscale = xscale = frame = 0;
    }

    void   set_viewport( double  _i_width,
                         double  _i_height,
                         double  _i_frame,
                         Box2d  & b_v )
    {
        width = _i_width;
        height = _i_height;
        frame = _i_frame;
        
        bb = b_v;
        xscale = (width - 2.0*frame) / (bb.xmax() - bb.xmin() );
        yscale = (height - 2.0*frame) / (bb.ymax() - bb.ymin() );
    }
    
    void   reset_viewport( Box2d   b_v )
    {
        bb = b_v;
        xscale = (width - 2.0*frame) / (bb.xmax() - bb.xmin() );
        //yscale = (height - 2.0*frame) / (bb.ymax() - bb.ymin() );

        double  yext;

        // yext = bb.ymax() - bb.ymin()
        yext = (height - 2.0*frame) / xscale;
        
        Point2d  cen = Point2d( b_v.xmin() + (bb.xmax() - bb.xmin())/2.0,
                                b_v.ymin() + (bb.ymax() - bb.ymin())/2.0 );

        double  y_min = cen.y() - yext/2.0;
        double  y_max = cen.y() + yext/2.0;
        
        bb = Box2d( bb.xmin(), y_min, bb.xmax(), y_max );
        
        yscale = xscale;
        
        //printf(  "xscale: %g, yscale: %g\n", xscale, yscale );
    }
    

    void  set_xscale( double  _xscale ) { xscale = _xscale; }
    void  set_yscale( double  _yscale ) { yscale = _yscale; }
    void  set_frame( int  _frame ) { frame = _frame; }

    double  scale() { return  xscale; }
    
    Point2d  map( double  x, double   y ) {
        double  rx, ry;
        
        rx =  ( x - bb.xmin() ) * xscale + frame;
        //ry = ( y - bb.ymin() ) * scale + frame;
        ry = ( bb.ymax() - y ) * yscale + frame;

        return  Point2d( rx, ry );
    }
    Point2d  map( const Point2d   & p ) {
        return  map( p.x(), p.y() );
    }
    double   get_in_height() const {
        return  bb.ymax() - bb.ymin() + 1;
    }
    double   get_in_width() const {
        return  bb.xmax() - bb.xmin() + 1;
    }
    
    int   get_out_width() const {
        return   2*frame + (int)(xscale * get_in_width());
    }
    int   get_out_height() const {
        return   2*frame + yscale * get_in_height();
    }
};

//////////////////////////////////////////////////////////////////////////
/// Utilities
//////////////////////////////////////////////////////////////////////////

void  bbox_print( const Box2d  & bb )
{
    printf( "BB: [%12f...%f] * [%f...%f]\n",
            bb.xmin(), bb.xmax(),
            bb.ymin(), bb.ymax() );
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

class  POI
{    
public:
    Point2d  loc;
    long  id;
    double  d;
    
    POI( long  _id, const double  & x, const double  & y ) {
        loc = Point2d( x, y );
        id = _id;
        d = 0;
    }
        
};
    

using  raw_edge_t = std::pair<long,long>;
typedef    cairo_t *cairo_t_ptr;
typedef    cairo_surface_t *cairo_surface_t_ptr;

class  GraphExt
{
private:
    vertices_id_map_t   v_id_map;
    int  n_vertices;
    std::vector<POI>  locs;
    std::vector<raw_edge_t>  raw_edges;
    Graph  g;
    Box2d  bbox;
    int  src;

    bool  f_png;
    int  page_drawn;
    
    std::vector<Vertex>  sp_p;
    std::vector<double>  sp_d;

    // Neck Cut values
    double totalDist;
    long edgeCountWavefront;
    double neckRatioMax = 0;
    double neckRatioDistMax = 0;

    void  draw_segment( cairo_t *cr, Mapper  & trans, Edge  e,
                        double  dist );
    void  draw_page( Mapper  & trans, cairo_t *cr, double  dist );
    void  draw_page_ext( Mapper      & trans,
                         cairo_surface_t_ptr  & surface,
                         cairo_t_ptr  & cr,
                         double  dist );

public:
    GraphExt() {
        n_vertices = 0;
        src = -1;
    }
    
    int  n() const { return  n_vertices; }

    Box2d  get_bounding_box();

    Box2d  get_active_bbox( double  dist );

    double  distance( int  a,  int  b )
    {
        const Point2d  & p_a = locs[ a ].loc;
        const Point2d  & p_b = locs[ b ].loc;

        return  sqrt(squared_distance( p_a, p_b ) );        
    }
    void  convert_raw_edges()
    {
        for  ( auto  p : raw_edges ) {
            int  a, b;

            //printf( "%ld - %ld\n", p.first, p.second );
            a = get_v_by_id(  p.first );
            b = get_v_by_id(  p.second );
            std::pair<Edge, bool>  ed = add_edge( a, b, g );
            Edge e = ed.first;
            
            // Set new edge weight to 0
            boost::put(boost::edge_weight_t(), g, e, 0 );
        }
    }

    void  update_edges_len()
    {
        EdgeIterator  ei, ei_end;

        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {            
            int  a, b;

            a = source(*ei, g);
            b = target(*ei, g);
            
            boost::put(boost::edge_weight_t(), g, *ei,
                       distance( a, b ) );
        }
    }
    

    
    void  insert_v_id( long  id, double  x, double  y )
    {
        static  double   y_min = 1000000000.0;
        static  double   y_max = -1;
        
        if ( v_id_map.find( id ) != v_id_map.end() )
            return;

        double  n_x = x * LONG_TO_METERS;
        double  n_y = y * LONG_TO_METERS;

        if  ( n_y < y_min )
            y_min = n_y;
        if  ( n_y > y_max )
            y_max = n_y;
        
        //printf( "Y: %12f...%12f\n", y_min, y_max );
        
        POI  p( id, n_x, n_y );
        
        locs.push_back( p ); //Point2d( x, y ) );
        v_id_map[ id ] = n_vertices++;        
    }

    int  get_v_by_id( long  id ) {
        auto  it = v_id_map.find( id );
        
        if ( it == v_id_map.end() ) {
            printf( "Error: Did not find %ld\n", id );
            return  -1;
        }
        return  it->second;
    }
    
    
    void  insert_raw_edge( long  ia, long ib ) {
        //printf( "Inserting edge...\n" );
        //printf( "ID = _%ld_\n", ia );
        //printf( "ID = _%ld_\n", ib );
        
        raw_edges.push_back( raw_edge_t( ia, ib ) );
        
        //add_edge( get_v_by_id( ia ),
        //           get_v_by_id( ib ),
        //  g );
    }
    
    void  shift_points_to_origin();

    void  draw( Mapper  & trans, const char  * out_pdf,
                double  dist );

    void    get_source( double  r_x, double  r_y );
    int    get_nearest_loc( const  Point2d   & p );
    void   compute_shortest_path();

    void   draw_dist_range( Mapper      & trans,
                            cairo_surface_t_ptr  & surface,
                            cairo_t_ptr  & cr,
                            double  d_start,
                            double  d_end,
                            double  d_delta );

};

int   GraphExt::get_nearest_loc( const  Point2d   & p )
{
    double  curr_d, d;
    int  curr;

    curr = 0;
    curr_d = squared_distance( p, locs[ 0 ].loc );
    for  ( int  ind = 1; ind < (int)locs.size(); ind++ ) {
        d = squared_distance( p, locs[ ind ].loc );
        if  ( d < curr_d ) {
            curr_d = d;
            curr = ind;
        }
    }

    return  curr;
}



void  GraphExt::get_source( double  r_x, double  r_y )
{
    bbox = get_bounding_box();

    Point2d  p( bbox.xmin() + (bbox.xmax() - bbox.xmin()) * r_x,
                bbox.ymin() + (bbox.ymax() - bbox.ymin()) * r_y );

    
    src = get_nearest_loc( p );
}


void  GraphExt::shift_points_to_origin()
{
    Box2d  bb;
    bb = get_bounding_box();

    
    for  ( int  ind = 0; ind < (int)locs.size(); ind++ ) {
        Point2d  p = locs[ ind ].loc;
        locs[ ind ].loc = Point2d( p.x() - bb.xmin(),
                                   p.y() - bb.ymin() );        
    }
        
}
        
void   bbox_bound_symmetric( Box2d  & b,
                             Point2d  & cen,
                             const Point2d  & nw )
{
    Box2d  bta( cen.x(), cen.y(),cen.x(), cen.y() );
    Box2d  btb( nw.x(), nw.y(), nw.x(), nw.y() );

    //Vector2d  v_nw =   nw - CGAL::ORIGIN;
    Point2d  ref = cen + (nw - cen) * (-1.0);
    
    Box2d  btc( ref.x() , ref.y(), ref.x(), ref.y() );
    
    b += bta;
    b += btb;
    b += btc;
}


Box2d  old_bbox;
bool  f_first = true;

Box2d  GraphExt::get_active_bbox( double  dist )
{
    Box2d  b;

    if  ( f_first )
        f_first = false;
    else {
        b = old_bbox;
    }
    Point2d  cen = locs[ src ].loc;

    /*
    for  ( int  ind = 0; ind < (int)locs.size(); ind++ ) {
        if  ( locs[ ind ].d <= dist )
            bbox_bound_symmetric( b, cen, locs[ ind ].loc );
    }

    
    Point2d  p( b.xmin() - frame, b.ymin() - frame );
    Point2d  q( b.xmin() + frame, b.ymin() + frame );

    bbox_bound_symmetric( b, cen, p );
    bbox_bound_symmetric( b, cen, q );
    */
     
    double  frame = 120;
    Point2d  pa( cen.x() - dist - frame, cen.y() - dist - frame );

    bbox_bound_symmetric( b, cen, pa );

    old_bbox = b;
    return  b;
}


Box2d  GraphExt::get_bounding_box()
{
    Box2d  b;

    for  ( auto loc : locs ) {
        Box2d  bt( loc.loc.x(), loc.loc.y(),
                   loc.loc.x(), loc.loc.y() );
        
        b += bt;

        //printf( "LOC: (%e, %e)\n",  loc.loc.x(), loc.loc.y() );
        //bbox_print( b );
    }
    //bbox_print( b );

    return  b;
}



class MyHandler : public osmium::handler::Handler {
public:
    GraphExt  * p_g;
    bool  f_reg_vertices;//, f_ins_edges;
    
    void way(const osmium::Way& way) {
        int  count = 0;
        long  prev = -1;
        
        //std::cout << "way " << way.id() << '\n';
        for (const auto& n : way.nodes()) {
            count++; 
            if  ( f_reg_vertices  &&  ( p_g != NULL ) )
                p_g->insert_v_id(  n.ref(), n.lon(), n.lat() );
                                   
            //std::cout << n.ref() << ": " << n.lon()
            //             << ", " << n.lat() << '\n';

            if  ( count > 1 ) 
                p_g->insert_raw_edge( prev, n.ref() );                
            
            prev = n.ref();
        }
    }
};


void  GraphExt::compute_shortest_path()
{
    
    //std::vector<vertex_descriptor> p(num_vertices(g));
    //std::vector<int> d(num_vertices(g));
    sp_p.resize( num_vertices(g), 0 );
    sp_d.resize( num_vertices(g), 0 );

    //property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
    //auto   pm =;

    std::vector<double> d(num_vertices(g));

    std::vector<double>  sp_d;
 
    dijkstra_shortest_paths( g, src,
      distance_map(boost::make_iterator_property_map(
            d.begin(), get(boost::vertex_index, g) )) );

    assert( locs.size() == d.size() );
    for  ( int  ind = 0; ind < (int)locs.size(); ind++ ) {
        locs[ ind ].d = d[ ind ];
    }
    
        
    /*
    VertexIterator v, vend;
    for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
        printf( "%d  dist: %g\n", (int)*v,
                boost::get(boost::vertex_distance_t(), g, *v ) );        
    }
    */
}


Point2d  convex_combination( const Point2d  & p_a,
                             const Point2d  & p_b,
                             double  delta )
{
    Point2d  p( p_a.x()  + delta * ( p_b.x() - p_a.x() ),
                p_a.y()  + delta * ( p_b.y() - p_a.y() ) );
    return  p;                 
}



void  GraphExt::draw_segment( cairo_t *cr,
                              Mapper  & trans,
                              Edge  e,
                              double  dist )
{
    cairo_set_line_cap  (cr, CAIRO_LINE_CAP_ROUND);

    Vertex v_a = source( e, g );
    Vertex v_b = target( e, g );

    double  len = distance( v_a, v_b );

    if ( locs[ v_a ].d > locs[ v_b ].d )
        std::swap( v_a, v_b );

    POI  & l_a( locs[ v_a ] );
    POI  & l_b( locs[ v_b ] );
        
    Point2d  pa, pb;
    
    pa = trans.map( l_a.loc );
    pb = trans.map( l_b.loc );
    
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width(cr, 2 );
    
    cairo_move_to(cr, pa.x(), pa.y() );
    cairo_line_to(cr, pb.x(), pb.y() );
    cairo_stroke(cr);

    /// Drawing the wter level...
    if ( l_a.d > dist )
        return;
    
    // double scale_param = 1.2;
    double scale_param = 1;

    // Edge is already traversed, color blue
    cairo_set_source_rgb( cr, 0, 0, 1 );
    cairo_set_line_width(cr, 8 );    
    if ( ( l_a.d + scale_param * len) <= dist ) {
        if (l_a.d<neckRatioDistMax && l_b.d>=neckRatioDistMax){
             cairo_set_source_rgb( cr, 0, 1, 0);
        }
        cairo_move_to(cr, pa.x(), pa.y() );
        cairo_line_to(cr, pb.x(), pb.y() );
        cairo_stroke(cr);
        cairo_set_source_rgb( cr, 0, 0, 1);
        totalDist += std::sqrt((pb.x()-pa.x())*(pb.x()-pa.x()) + (pb.y()-pa.y())*(pb.y()-pa.y()));
        return;
    }

    /// We have to draw the water on the segment from both sides. Yuk

    /// We start with the case that (we imagine that v_a is inside,
    /// but v_b is outside.

    // We are traversing this edge from v_a.
    double  delta = ( dist - l_a.d ) / len; // between zero and one.
    if  ( delta > 1.0 )
        delta = 1.0;
    
    Point2d  mid = convex_combination( l_a.loc, l_b.loc, delta );
    Point2d  pmid = trans.map( mid );
    
    cairo_move_to(cr, pa.x(), pa.y() );
    cairo_line_to(cr, pmid.x(), pmid.y() );
    cairo_stroke(cr);
    double drawn_len = std::sqrt((pmid.x()-pa.x())*(pmid.x()-pa.x()) + (pmid.y()-pa.y())*(pmid.y()-pa.y()));

    if ( l_b.d > dist ){
        edgeCountWavefront+=1;
        totalDist += drawn_len;
        return;
    }

    /// We have to draw the water from v_b...
    delta = ( dist - l_b.d ) / len; // between zero and one.

    mid = convex_combination( l_b.loc, l_a.loc, delta );
    pmid = trans.map( mid );
    cairo_move_to(cr, pb.x(), pb.y() );
    cairo_line_to(cr, pmid.x(), pmid.y() );
    cairo_stroke(cr);
    drawn_len +=  std::sqrt((pmid.x()-pb.x())*(pmid.x()-pb.x()) + (pmid.y()-pb.y())*(pmid.y()-pb.y()));
    totalDist += drawn_len;
    edgeCountWavefront+=2;
}


void  GraphExt::draw_page( Mapper  & trans,
                           cairo_t *cr,
                           double  dist )
{

    cairo_set_line_cap  (cr, CAIRO_LINE_CAP_ROUND);

    /* Set surface to opaque color (r, g, b) */
    cairo_set_source_rgb (cr, 1, 1, 1 );
    cairo_paint (cr);
    //cairo_stroke(cr);
    page_drawn++;   // count the pages already drawn
    
    // Zoom out as distance increases
    Box2d  n_bb = get_active_bbox( dist * 1.5 );
    
    trans.reset_viewport( n_bb );

    
    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_set_line_width(cr, 0.2 );    
    //////////////////////////////////////////////////////////

    // Draw every edge in the graph.
    edgeCountWavefront = 0;
    totalDist = 0;
    EdgeIterator  ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        cairo_set_line_cap  (cr, CAIRO_LINE_CAP_ROUND);
        draw_segment( cr, trans,  *ei, dist );
        
    }

    // Draw source point
    if  ( src >= 0 ) {
        Point2d   m = trans.map( locs[ src ].loc );
        
        cairo_set_source_rgb(cr, 1, 0, 0);
        ///cairo_arc(cr, m.x(), m.y(), 5.0 * trans.scale(), 0, 2 * M_PI);
        cairo_arc(cr, m.x(), m.y(), 10.0, 0, 2 * M_PI);
        cairo_fill(cr);
        cairo_stroke(cr);            
    }

    /// Write the distance
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_BOLD);
    
    cairo_set_font_size(cr, 24 );

    cairo_set_source_rgb(cr,0,0,0);
    cairo_set_line_width(cr, 2);    
   
    char  buf[ 256 ];
    sprintf( buf, "Distance: %g  ", dist );


    cairo_move_to(cr, 20, 34);
    cairo_show_text(cr, buf );
    cairo_stroke(cr);

    // Print accumulated distance
    sprintf( buf, "Total Distance: %g  ", totalDist );


    cairo_move_to(cr, 20, 68);
    cairo_show_text(cr, buf );
    cairo_stroke(cr); 

    // Print Edges in wavefront
    sprintf( buf, "Current Wavefront: %ld  ", edgeCountWavefront );


    cairo_move_to(cr, 20, 34*3);
    cairo_show_text(cr, buf );
    cairo_stroke(cr);

    // Compute Ratio, total dist / wavefront count **2 
    double neck_ratio =   totalDist / (double)(edgeCountWavefront*edgeCountWavefront);
    if ( std::max(neck_ratio, neckRatioMax)==neck_ratio){
        neckRatioMax = neck_ratio;
        neckRatioDistMax = dist;
    }

    // Neck Cut ratio
    sprintf( buf, "CNCR: %g; MNCR: %g, @ Dist: %g", neck_ratio, neckRatioMax, neckRatioDistMax);


    cairo_move_to(cr, 20, 34*4);
    cairo_show_text(cr, buf );
    cairo_stroke(cr);

    
    cairo_show_page( cr );
}

typedef    cairo_t *cairo_t_ptr;

void  GraphExt::draw_page_ext( Mapper      & trans,
                               cairo_surface_t_ptr & surface,
                               cairo_t_ptr  & cr,
                               double  dist )
{
    char   buf[ 1024 ];
    
    // Set up cairo surface for images
    if  ( f_png ) {
        sprintf( buf, "frames/%06d.png", page_drawn );
        printf( "  %s\n", buf );
        surface = cairo_image_surface_create
            ( CAIRO_FORMAT_RGB24, trans.get_out_width(),
              trans.get_out_height() );
        cr = cairo_create(surface);
    }
    
    // Draw image
    draw_page( trans, cr, dist );
    
    // Write to png and clean up surface.
    if  ( f_png ) {
        cairo_surface_write_to_png( surface, buf );
        cairo_destroy (cr);
        cairo_surface_destroy (surface);
    }
}

// This function draws a frame for each distance value in range, incremented by delta
void  GraphExt::draw_dist_range(
                                Mapper      & trans,
                                cairo_surface_t_ptr  & surface,
                                cairo_t_ptr  & cr,
                                double  d_start,
                                double  d_end,
                                double  d_delta )
{
    double  d;

    d = d_start;

    // Draw a page for each dist in our range, increasing it by delta each iteration.
    do { 
        draw_page_ext( trans, surface, cr, d );
        d += d_delta;
    } while  ( d < d_end );    
}

// This function draws the graph + continous dijkstra as as series of pngs, or onto a page on the output pdf.                         

void  GraphExt::draw( Mapper  & trans, const char  * out_pdf,
                      double  dist )
{
    // Should we output png frames
    f_png = true;

    //////////////////////////////////////////////////////////
    // Cairo setup
    //////////////////////////////////////////////////////////
    cairo_surface_t *surface;
    cairo_t *cr;

    // Frame or page counter
    page_drawn = 0;

    if  ( ! f_png ) {
        //Create PDF surface from cairo
        surface = (cairo_surface_t *)cairo_pdf_surface_create
        ( out_pdf, trans.get_out_width(),
          trans.get_out_height() );
        cr = cairo_create(surface);
    } else {
        // If we're not using pdf builder, NULL them
        surface = NULL;
        cr = NULL;
    }

    // Draw the series of distances with increasing rates as we grow further from source.
    draw_dist_range( trans, surface, cr, 0, 200, 1 );
    draw_dist_range( trans, surface, cr, 200, 400, 2 );
    draw_dist_range( trans, surface, cr, 400, 800, 2 );
    draw_dist_range( trans, surface, cr, 800, 1600, 4 );
    draw_dist_range( trans, surface, cr, 1600, 3200, 8 );
    draw_dist_range( trans, surface, cr, 3200, 6400, 32 );

    // Clean up cairo surfaces if outputing to pdf.
    if  ( ! f_png ) {
        cairo_destroy (cr);
        cairo_surface_destroy (surface);
    }    
}



void  read_map( GraphExt  & ge, const char  * filename )
{
    auto otypes = osmium::osm_entity_bits::node | osmium::osm_entity_bits::way;
    (void)otypes;
    osmium::io::Reader reader{ filename };

    namespace map = osmium::index::map;
    using index_type = map::SparseMemArray<osmium::unsigned_object_id_type,
                                           osmium::Location>;
    using location_handler_type
        = osmium::handler::NodeLocationsForWays<index_type>;

    index_type index;
    location_handler_type location_handler{index};

    MyHandler handler;
    handler.p_g = &ge;
    handler.f_reg_vertices = true;
    //handler.f_ins_edges = false;
        
    osmium::apply(reader, location_handler, handler);
    //osmium::apply(reader, location_handler, handler);


    reader.close();

    printf( "Converting raw edges...\n" );
    
    ge.convert_raw_edges();
        
    printf( "Number of vertices: %d\n", ge.n() );
}




int main(int argc, char* argv[])
{
    GraphExt  ge;
    double image_width = 1920;
    double image_height = 1080;
    double image_frame = 10;


    // Choose a default map, or one supplied by cmd line.
    if (argc==1){
        //read_map( ge, "data/highways.osm" );
        read_map( ge, "data/siebel.osm" );
        // read_map( ge, "data/old_jerusalem_r.osm" );
    } else{
        read_map( ge, argv[1]);
    }
    
    
    ge.shift_points_to_origin();
    ge.update_edges_len();

    
    ////////////////////////////////////////////////////////////////    
    Box2d  bb;
    bb = ge.get_bounding_box();
        
    bbox_print( bb );
    Mapper  trans( bb );

    trans.set_viewport( image_width, image_height, image_frame, bb );
    
    //double xscale, yscale;


    printf( "width  = %d\n", trans.get_out_width() );
    printf( "height = %d\n", trans.get_out_height() );

    //ge.get_source( 0.75, 0.5 );
    ge.get_source( 0.45, 0.53 );
    ge.compute_shortest_path();
    
    ge.draw( trans, "map.pdf", 200 );
    
    return  0;
}
