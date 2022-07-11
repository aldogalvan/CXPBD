//
// Created by aldo on 6/10/22.
//

#ifndef CXPBD_CEDGEARRAY_H
#define CXPBD_CEDGEARRAY_H

#include "chai3d.h"
namespace chai3d{

//------------------------------------------------------------------------------
class cEdgeArray;
typedef std::shared_ptr<cEdgeArray> cEdgeArrayPtr;
//------------------------------------------------------------------------------

//==============================================================================
class cEdgeArray : public cGenericArray
{
    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //--------------------------------------------------------------------------
    /*!
        Constructor of cEdgeArray.

        \param  a_vertexArray  Array of vertices used to describe the edge.
    */
    //--------------------------------------------------------------------------
    cEdgeArray(cVertexArrayPtr a_vertexArray) : cGenericArray(a_vertexArray)
    {
        // clear all triangles
        clear();

        // store pointer to vertex array
        m_vertices = a_vertexArray;

        // initialize OpenGL buffer ID
        m_elementBuffer         = 0;
    }


    //--------------------------------------------------------------------------
    /*!
        Destructor of cEdgeArray.
    */
    //--------------------------------------------------------------------------
    ~cEdgeArray(){}


    //! Clear all edges from array.
    void clear()
    {
        m_allocated.clear();
        m_indices.clear();
        m_freeElements.clear();
        m_flagMarkForUpdate     = true;
        m_flagMarkForResize     = true;
    }


    //! Shared cTriangleArray allocator.
    static cEdgeArrayPtr create(cVertexArrayPtr a_vertexArray) { return (std::make_shared<cEdgeArray>(a_vertexArray)); }


    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:

    //! This method returns the number of vertices that compose a edge.
    virtual unsigned int getNumVerticesPerElement() { return (2); }

    //! This method copies all edge data and returns a new edge array.
    cEdgeArrayPtr copy();

    //! This method removes all non allocated triangles.
    void compress();


    //--------------------------------------------------------------------------
    /*!
        This method creates a new edge.

        \param  a_vertexIndex0  index of vertex 0.
        \param  a_vertexIndex1  index of vertex 1.


        \return Index number of new edge.
    */
    //--------------------------------------------------------------------------
    int newEdge(const unsigned int a_vertexIndex0,
                       const unsigned int a_vertexIndex1)
    {
        int index;

        // check if there is an available slot on the free tetrahedra list
        if (m_freeElements.size() > 0)
        {
            index = m_freeElements.front();
            m_allocated[index] = true;
            m_freeElements.pop_front();
            setVertices(index, a_vertexIndex0, a_vertexIndex1);
        }
        else
        {
            // store new indices
            m_indices.push_back(a_vertexIndex0);
            m_indices.push_back(a_vertexIndex1);
            m_allocated.push_back(true);

            m_flagMarkForResize = true;

            // store index of new edge
            index = (int)(m_allocated.size())-1;
        }

        // mark for update
        m_flagMarkForUpdate = true;

        // return index to new edge
        return (index);
    }


    //--------------------------------------------------------------------------
    /*!
        This method deallocates a selected edge from the array. The three
        vertices of the edge are set to zero, and the edge is added to
        the free list for future allocation.\n

        \param  a_edgeIndex  Index of selected edge.
    */
    //--------------------------------------------------------------------------
    void removeEdge(const unsigned int a_edgeIndex)
    {
        // sanity check
        if (m_allocated[a_edgeIndex])
        {
            // disable triangle
            m_allocated[a_edgeIndex] = false;

            // reset indices to vertices
            setVertices(a_edgeIndex, 0, 0);

            // add triangle to free list
            m_freeElements.push_back(a_edgeIndex);

            // mark for update
            m_flagMarkForUpdate = true;
        }
    }


    //--------------------------------------------------------------------------
    /*!
        This method set the vertices of a selected edge by passing their
        index numbers.

        \param  a_tetrahedronIndex  Index of selected edge.
        \param  a_vertexIndex0   Index of vertex 0.
        \param  a_vertexIndex1   Index of vertex 1.


    */
    //--------------------------------------------------------------------------
    inline void setVertices(const unsigned int a_tetrahedronIndex,
                            const unsigned int a_vertexIndex0,
                            const unsigned int a_vertexIndex1)
    {
        m_indices[2*a_tetrahedronIndex+0] = a_vertexIndex0;
        m_indices[2*a_tetrahedronIndex+1] = a_vertexIndex1;

        // mark for update
        m_flagMarkForUpdate = true;
    }


    //--------------------------------------------------------------------------
    /*!
        This method sets vertex 0 of a selected edge by passing their
        index numbers.

        \param  a_triangleIndex  Index of edge.
        \param  a_vertexIndex0   Index of vertex 0.
    */
    //--------------------------------------------------------------------------
    inline void setVertexIndex0(const unsigned int a_edgeIndex, const unsigned int a_vertexIndex0)
    {
        m_indices[2 * a_edgeIndex + 0] = a_vertexIndex0;

        // mark for update
        m_flagMarkForUpdate = true;
    };


    //--------------------------------------------------------------------------
    /*!
        This method returns the index of vertex 0 of a selected edge.

        \param  a_edgeIndex  Index of edge.

        \return Index to vertex 0.
    */
    //--------------------------------------------------------------------------
    inline unsigned int getVertexIndex0(const unsigned int a_edgeIndex) const
    {
        return (m_indices[2*a_edgeIndex+0]);
    };


    //--------------------------------------------------------------------------
    /*!
        This method sets vertex 1 of a selected edge by passing their
        index numbers.

        \param  a_edgeIndex  Index of edge.
        \param  a_vertexIndex1   Index of vertex 1.
    */
    //--------------------------------------------------------------------------
    inline void setVertexIndex1(const unsigned int a_edgeIndex, const unsigned int a_vertexIndex1)
    {
        m_indices[2 * a_edgeIndex + 1] = a_vertexIndex1;

        // mark for update
        m_flagMarkForUpdate = true;
    };


    //--------------------------------------------------------------------------
    /*!
        This method returns the index of vertex 1 of a selected edge.

        \param  a_edgeIndex  Index of selected edge.

        \return Index to vertex 1.
    */
    //--------------------------------------------------------------------------
    inline unsigned int getVertexIndex1(const unsigned int a_edgeIndex) const
    {
        return (m_indices[2*a_edgeIndex+1]);
    };



    //--------------------------------------------------------------------------
    /*!
         This method returns the index number of a selected vertex for a
         selected edge. In the case of a edge a_vertexNumber takes
         either 0, or 2.

        \param  a_elementIndex  Index number of selected edge.
        \param  a_vertexNumber  Vertex number.

        \return Index number of the requested vertex.
    */
    //--------------------------------------------------------------------------
    virtual unsigned int getVertexIndex(const unsigned int a_elementIndex,
                                        const unsigned int a_vertexNumber) const
    {
        switch (a_vertexNumber)
        {
            case 0: return (m_indices[4*a_elementIndex+0]);
            case 1: return (m_indices[4*a_elementIndex+1]);
            default: return (0);
        }
    }


    //--------------------------------------------------------------------------
    /*!
        This method renders the OpenGL vertex buffer object.
    */
    //--------------------------------------------------------------------------
    inline void renderInitialize()
    {
#ifdef C_USE_OPENGL
        unsigned int numtetrahedra = getNumElements();

        // allocate object and buffers if needed
        if (m_elementBuffer == 0)
        {
            glGenBuffers(1, &m_elementBuffer);
        }

        // bind buffer
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_elementBuffer);

        // update buffer size if needed
        if (m_flagMarkForResize)
        {
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 2 * numedges * sizeof(unsigned int), &(m_indices[0]), GL_STATIC_DRAW);
            m_flagMarkForResize = true;
        }

        // update data if needed
        if (m_flagMarkForUpdate)
        {
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_elementBuffer);
            glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, 2 * numedges * sizeof(unsigned int), &(m_indices[0]));

            m_flagMarkForUpdate = false;
        }

        // initialize rendering of vertices
        m_vertices->renderInitialize();

        // render object
        glEnableClientState(GL_VERTEX_ARRAY);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_elementBuffer);
        glDrawElements(GL_EDGES, 2 * numedges, GL_UNSIGNED_INT, (void*)0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
#endif
    }

    //--------------------------------------------------------------------------
    /*!
        This method finalizes rendering the OpenGL vertex buffer object.
    */
    //--------------------------------------------------------------------------
    inline void renderFinalize()
    {
#ifdef C_USE_OPENGL
        // finalize rendering
        m_vertices->renderFinalize();
        glDisableClientState(GL_VERTEX_ARRAY);
#endif
    }



    //--------------------------------------------------------------------------
    /*!
        This method computes and returns the length for a selected edge.

        \return length for selected edge.
    */
    //--------------------------------------------------------------------------
    inline double computeLength(const unsigned int a_edgeIndex)
    {
        // A = 0.5 * | u x v |
        cVector3d a = m_vertices->getLocalPos(getVertexIndex1(a_edgeIndex));
        cVector3d b = m_vertices->getLocalPos(getVertexIndex2(a_edgeIndex));


        return (cSub(a,b).length());
    }




//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

#endif //CXPBD_CEDGEARRAY_H
