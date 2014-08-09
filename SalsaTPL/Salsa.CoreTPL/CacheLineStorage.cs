using System.Runtime.InteropServices;

namespace Salsa.Core
{
    // C#: Note works for value types only
    // It may seem strange that this code actually allocates enough space for two cache lines' 
    // worth of data instead of just one. That's because, on .NET, you can't specify the alignment 
    // of data beyond some inherent 4-byte and 8-byte alignment guarantees, which aren't big enough 
    // for our purposes. Even if you could specify a starting alignment, the compacting garbage 
    // collector is likely to move your object and thus change its alignment dynamically. Without 
    // alignment to guarantee the starting address of the data, the only way to deal with this is 
    // to allocate enough space both before and after data to ensure that no other objects can 
    // share the cache line.
    [StructLayout(LayoutKind.Explicit, Size = 128)]
    public struct CacheAwareStorage<T> where T : struct
    {
        [FieldOffset(64)] public T Value;
    }

    public class CacheAwareArray<T> where T : struct
    {
        private readonly CacheAwareStorage<T>[] _items;

        public CacheAwareArray(int length)
        {
            _items = new CacheAwareStorage<T>[length + 1];
        }

        public T this[int index]
        {
            get { return _items[index + 1].Value; }
            set { _items[index + 1].Value = value; }
        }

        public int Length
        {
            get { return _items.Length - 1; }
        }
    }
}