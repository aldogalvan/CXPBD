//
// Created by aldo on 6/10/22.
//

#ifndef CXPBD_CTETRAHEDRAARRAY_H
#define CXPBD_CTETRAHEDRAARRAY_H
#include "chai3d.h"
namespace chai3d
{

//------------------------------------------------------------------------------
class cTetrahedraArray;
typedef std::shared_ptr<cTetrahedraArray> cTetrahedraArrayPtr;
//------------------------------------------------------------------------------

//==============================================================================
class cTetrahedraArray : public cGenericArray
{
    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //--------------------------------------------------------------------------
    /*!
        Constructor of cTriangleArray.

        \param  a_vertexArray  Array of vertices used to describe the triangles.
    */
    //--------------------------------------------------------------------------
    cTetrahedraArray(cVertexArrayPtr a_vertexArray) : cGenericArray(a_vertexArray)
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
        Destructor of cTriangleArray.
    */
    //--------------------------------------------------------------------------
    ~cTetrahedraArray(){}


    //! Clear all triangles from array.
    void clear()
    {
        m_allocated.clear();
        m_indices.clear();
        m_freeElements.clear();
        m_flagMarkForUpdate     = true;
        m_flagMarkForResize     = true;
    }


    //! Shared cTriangleArray allocator.
    static cTetrahedraArrayPtr create(cVertexArrayPtr a_vertexArray) { return (std::make_shared<cTetrahedraArray>(a_vertexArray)); }


    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:

    //! This method returns the number of vertices that compose a triangle.
    virtual unsigned int getNumVerticesPerElement() { return (4); }

    //! This method copies all triangle data and returns a new triangle array.
    cTetrahedraArrayPtr copy();

    //! This method removes all non allocated triangles.
    void compress();


    //--------------------------------------------------------------------------
    /*!
        This method creates a new triangle.

        \param  a_vertexIndex0  index of vertex 0.
        \param  a_vertexIndex1  index of vertex 1.
        \param  a_vertexIndex2  index of vertex 2.
        \param  a_vertexIndex3  index of vertex 3

        \return Index number of new tetrahedra.
    */
    //--------------------------------------------------------------------------
    int newTetrahedron(const unsigned int a_vertexIndex0,
                       const unsigned int a_vertexIndex1,
                       const unsigned int a_vertexIndex2,
                       const unsigned int a_vertexIndex3)
    {
        int index;

        // check if there is an available slot on the free tetrahedra list
        if (m_freeElements.size() > 0)
        {
            index = m_freeElements.front();
            m_allocated[index] = true;
            m_freeElements.pop_front();
            setVertices(index, a_vertexIndex0, a_vertexIndex1,
                        a_vertexIndex2, a_vertexIndex3);
        }
        else
        {
            // store new indices
            m_indices.push_back(a_vertexIndex0);
            m_indices.push_back(a_vertexIndex1);
            m_indices.push_back(a_vertexIndex2);
            m_indices.push_back(a_vertexIndex3)
            m_allocated.push_back(true);

            m_flagMarkForResize = true;

            // store index of new triangle
            index = (int)(m_allocated.size())-1;
        }

        // mark for update
        m_flagMarkForUpdate = true;

        // return index to new triangle
        return (index);
    }


    //--------------------------------------------------------------------------
    /*!
        This method deallocates a selected triangle from the array. The three
        vertices of the triangle are set to zero, and the triangle is added to
        the free list for future allocation.\n

        \param  a_triangleIndex  Index of selected triangle.
    */
    //--------------------------------------------------------------------------
    void removeTetrahedron(const unsigned int a_tetrahedronIndex)
    {
        // sanity check
        if (m_allocated[a_tetrahedronIndex])
        {
            // disable triangle
            m_allocated[a_tetrahedronIndex] = false;

            // reset indices to vertices
            setVertices(a_tetrahedronIndex, 0, 0, 0, 0);

            // add triangle to free list
            m_freeElements.push_back(a_tetrahedronIndex);

            // mark for update
            m_flagMarkForUpdate = true;
        }
    }


    //--------------------------------------------------------------------------
    /*!
        This method set the vertices of a selected tetrahedron by passing their
        index numbers.

        \param  a_tetrahedronIndex  Index of selected tetrahedron.
        \param  a_vertexIndex0   Index of vertex 0.
        \param  a_vertexIndex1   Index of vertex 1.
        \param  a_vertexIndex2   Index of vertex 2.
        \param  a_vertexIndex3   Index of vertex 3.

    */
    //--------------------------------------------------------------------------
    inline void setVertices(const unsigned int a_tetrahedronIndex,
                            const unsigned int a_vertexIndex0,
                            const unsigned int a_vertexIndex1,
                            const unsigned int a_vertexIndex2,
                            const unsigned int a_vertexIndex3)
    {
        m_indices[4*a_tetrahedronIndex+0] = a_vertexIndex0;
        m_indices[4*a_tetrahedronIndex+1] = a_vertexIndex1;
        m_indices[4*a_tetrahedronIndex+2] = a_vertexIndex2;
        m_indices[4*a_tetrahedronIndex+3] = a_vertexIndex2;

        // mark for update
        m_flagMarkForUpdate = true;
    }


    //--------------------------------------------------------------------------
    /*!
        This method sets vertex 0 of a selected tetrahedron by passing their
        index numbers.

        \param  a_triangleIndex  Index of tetrahederon.
        \param  a_vertexIndex0   Index of vertex 0.
    */
    //--------------------------------------------------------------------------
    inline void setVertexIndex0(const unsigned int a_tetrahedronIndex, const unsigned int a_vertexIndex0)
    {
        m_indices[4 * a_tetrahedronIndex + 0] = a_vertexIndex0;

        // mark for update
        m_flagMarkForUpdate = true;
    };


    //--------------------------------------------------------------------------
    /*!
        This method returns the index of vertex 0 of a selected tetrahedron.

        \param  a_tetrahedronIndex  Index of tetrahedron.

        \return Index to vertex 0.
    */
    //--------------------------------------------------------------------------
    inline unsigned int getVertexIndex0(const unsigned int a_tetrahedronIndex) const
    {
        return (m_indices[4*a_tetrahedronIndex+0]);
    };


    //--------------------------------------------------------------------------
    /*!
        This method sets vertex 1 of a selected tetrahedron by passing their
        index numbers.

        \param  a_tetrahedronIndex  Index of tetrahedron.
        \param  a_vertexIndex1   Index of vertex 1.
    */
    //--------------------------------------------------------------------------
    inline void setVertexIndex1(const unsigned int a_tetrahedronIndex, const unsigned int a_vertexIndex1)
    {
        m_indices[4 * a_tetrahedronIndex + 1] = a_vertexIndex1;

        // mark for update
        m_flagMarkForUpdate = true;
    };


    //--------------------------------------------------------------------------
    /*!
        This method returns the index of vertex 1 of a selected triangle.

        \param  a_triangleIndex  Index of selected triangle.

        \return Index to vertex 1.
    */
    //--------------------------------------------------------------------------
    inline unsigned int getVertexIndex1(const unsigned int a_tetrahedronIndex) const
    {
        return (m_indices[4*a_tetrahedronIndex+1]);
    };


    //--------------------------------------------------------------------------
    /*!
        This method sets vertex 2 of a selected triangle by passing their
        index numbers.

        \param  a_triangleIndex  Index of triangle.
        \param  a_vertexIndex2   Index of vertex 2.
    */
    //--------------------------------------------------------------------------
    inline void setVertexIndex2(const unsigned int a_tetrahedronIndex, const unsigned int a_vertexIndex2)
    {
        m_indices[4 * a_tetrahedronIndex + 2] = a_vertexIndex2;

        // mark for update
        m_flagMarkForUpdate = true;
    };


    //--------------------------------------------------------------------------
    /*!
        This method returns the index of vertex 2 of a selected triangle.

        \param  a_triangleIndex  Index of triangle.

        \return Index to vertex 2.
    */
    //--------------------------------------------------------------------------
    inline unsigned int getVertexIndex2(const unsigned int a_tetrahedronIndex) const
    {
        return (m_indices[4*a_tetrahedronIndex+2]);
    };

    //--------------------------------------------------------------------------
    /*!
        This method sets vertex 2 of a selected triangle by passing their
        index numbers.

        \param  a_triangleIndex  Index of triangle.
        \param  a_vertexIndex2   Index of vertex 2.
    */
    //--------------------------------------------------------------------------
    inline void setVertexIndex3(const unsigned int a_tetrahedronIndex, const unsigned int a_vertexIndex3)
    {
        m_indices[4 * a_tetrahedronIndex + 3] = a_vertexIndex3;

        // mark for update
        m_flagMarkForUpdate = true;
    };


    //--------------------------------------------------------------------------
    /*!
        This method returns the index of vertex 2 of a selected triangle.

        \param  a_triangleIndex  Index of triangle.

        \return Index to vertex 2.
    */
    //--------------------------------------------------------------------------
    inline unsigned int getVertexIndex3(const unsigned int a_tetrahedronIndex) const
    {
        return (m_indices[4*a_tetrahedronIndex+3]);
    };


    //--------------------------------------------------------------------------
    /*!
         This method returns the index number of a selected vertex for a
         selected triangle. In the case of a triangle a_vertexNumber takes
         either 0, 1 or 2.

        \param  a_elementIndex  Index number of selected triangle.
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
            case 2: return (m_indices[4*a_elementIndex+2]);
            case 3: return (m_indices[4*a_elementIndex+3]);
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
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 4 * numtetrahedra * sizeof(unsigned int), &(m_indices[0]), GL_STATIC_DRAW);
            m_flagMarkForResize = true;
        }

        // update data if needed
        if (m_flagMarkForUpdate)
        {
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_elementBuffer);
            glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, 4 * numtetrahedra * sizeof(unsigned int), &(m_indices[0]));

            m_flagMarkForUpdate = false;
        }

        // initialize rendering of vertices
        m_vertices->renderInitialize();

        // render object
        glEnableClientState(GL_VERTEX_ARRAY);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_elementBuffer);
        glDrawElements(GL_TRIANGLES, 4 * numtetrahedra, GL_UNSIGNED_INT, (void*)0);
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
        This method computes and returns the volume for a selected tetrahedron.

        \return volume for selected tetrahedron.
    */
    //--------------------------------------------------------------------------
    inline double computeVolume(const unsigned int a_tetrahedronIndex)
    {
        // A = 0.5 * | u x v |
        cVector3d u = m_vertices->getLocalPos(getVertexIndex1(a_tetrahedronIndex)) - m_vertices->getLocalPos(getVertexIndex0(a_tetrahedronIndex));
        cVector3d v = m_vertices->getLocalPos(getVertexIndex2(a_tetrahedronIndex)) - m_vertices->getLocalPos(getVertexIndex0(a_tetrahedronIndex));
        cVector3d z = m_vertices->getLocalPos(getVertexIndex3(a_tetrahedronIndex)) - m_vertices->getLocalPos(getVertexIndex0(a_tetrahedronIndex));

        return (1./6. * (cCross(u,v).cDot(z)));
    }




//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

#endif //CXPBD_CTETRAHEDRAARRAY_H