using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;

namespace Salsa.Core.Blas
{
    public class Vector<T> : IEnumerable<T>
    {
        protected T[] _elements;
        protected int _length;

        public Vector(T[] values)
        {
            _length = values.Length;
            _elements = values;
        }

        public Vector(int length)
        {
            _length = length;
            _elements = new T[length];
        }

        public T this[int index]
        {
            get { return _elements[index]; }
            set { _elements[index] = value; }
        }

        public int Length
        {
            get { return _length; }
        }

        #region IEnumerable<T> Members

        public IEnumerator<T> GetEnumerator()
        {
            for (int i = 0; i < _length; i++)
            {
                yield return _elements[i];
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        #endregion

        public virtual void SetAll(T value)
        {
            for (int i = 0; i < _length; i++)
            {
                _elements[i] = value;
            }
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            for (int i = 0; i < _elements.Length; i++)
            {
                if (i != 0)
                {
                    sb.Append(',');
                }

                sb.Append(_elements[i]);
            }

            return sb.ToString();
        }

        public static implicit operator T[](Vector<T> v)
        {
            return v._elements;
        }

        public static implicit operator Vector<T>(T[] v)
        {
            return new Vector<T>(v);
        }

        internal static bool CheckDimensions(Vector<T> left, Vector<T> right, bool throwException)
        {
            if (left._length != right._length)
            {
                if (throwException)
                {
                    throw new IndexOutOfRangeException();
                }

                return false;
            }

            return true;
        }
    }
}