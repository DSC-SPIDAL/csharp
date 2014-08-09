namespace Salsa.Core
{
    public delegate T OnDemandCompute<T>();

    public class OnDemandComputation<T>
    {
        private OnDemandCompute<T> _compute;
        private T _result;
        private bool _shouldCompute;

        public OnDemandComputation(OnDemandCompute<T> compute)
        {
            _compute = compute;
            _shouldCompute = true;
        }

        public bool ShouldCompute
        {
            get { return _shouldCompute; }
        }

        public T Compute()
        {
            if (_shouldCompute)
            {
                _result = Compute();
                _shouldCompute = false;
            }

            return _result;
        }

        public void Reset()
        {
            _result = default(T);
            _shouldCompute = true;
        }
    }
}